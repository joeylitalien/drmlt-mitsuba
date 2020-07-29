/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/bidir/util.h>
#include <mitsuba/core/plugin.h>
#include "drmlt_proc.h"
#include "drmlt_sampler.h"
#include "../blockthread.h"
#include <mutex>

MTS_NAMESPACE_BEGIN

/*!\plugin{drmlt}{Delayed Rejection Metropolis Light Transport}
* \order{9}
 * \parameters{
 *      \parameter{technique}{PathSampler::ETechnique}{
 *          DRMLT works in conjunction with another rendering
 *          technique that is endowed with Markov Chain-based sample generation.
 *          Three choices are available:
 *          \begin{itemize}
 *              \item \code{bdpt}: Operate on top of a fully-fleged bidirectional
 *                  path tracer with multiple importance sampling.
 *              \item \code{path}: Rely on a unidirectional
 *                  volumetric path tracer (i.e. \pluginref{volpath}) 
 *              \item \code{mmlt}: Rely on multiplexed Metropolis light transport.
 *                  Note that the stragtegy is kept fixed in this application for simplicity
 *              \vspace{-4mm}
 *          \end{itemize}
 *      }
 *      \parameter{type}{DRMLTConfiguration::ETechnique}{
 *          Three different DRMLT algorithm are available
 *          \begin{itemize}
 *              \item \code{green}: Use Green & Mira (2001) generalized delayed rejection 
 *                  framework. 
 *              \item \code{mira}: Use Tierney & Mira (1999) vanilla delayed rejection 
 *                  framework. 
 *              \item \code{orbital}: Use Tierney & Mira (1999) framework coupled with an
 *                  orbital mutation at the second stage. 
 *              \vspace{-4mm}
 *          \end{itemize}
 *      }
 *      \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *          in the generated output image (where \code{-1} corresponds to $\infty$).
 *	        A value of \code{1} will only render directly visible light sources.
 *	        \code{2} will lead to single-bounce (direct-only) illumination,
 *	        and so on. \default{\code{-1}}
 *	    }
 *      \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	        which the implementation will start to use the ``russian roulette''
 *          path termination criterion. \default{\code{5}}
 *	    }
 *      \parameter{directSampling}{\Boolean}{Enable direct sampling strategies? This is a generalization
 *          of direct illumination sampling that works with both emitters and sensors. Usually a good idea.
 *          \default{use direct sampling, i.e. \code{true}}
 *      }
 *	    \parameter{directSamples}{\Integer}{
 *	        By default, this plugin renders the direct illumination component
 *	        separately using an optimized direct illumination sampling strategy
 *	        that uses low-discrepancy number sequences for superior performance
 *	        (in other words, it is \emph{not} rendered by PSSMLT). This
 *	        parameter specifies the number of samples allocated to that method. To
 *	        force PSSMLT to be responsible for the direct illumination
 *	        component as well, set this parameter to \code{-1}. \default{16}
 *	    }
 *	    \parameter{luminanceSamples}{\Integer}{
 *	        MLT-type algorithms create output images that are only
 *	        \emph{relative}. The algorithm can e.g. determine that a certain pixel
 *	        is approximately twice as bright as another one, but the absolute
 *	        scale is unknown. To recover it, this plugin computes
 *	        the average luminance arriving at the sensor by generating a
 *	        number of samples. \default{\code{100000} samples}
 *      }
 *      \parameter{twoStage}{\Boolean}{Use two-stage MLT?
 *          See below for details. \default{{\footnotesize\code{false}}}
 *      }
 *	    \parameter{pLarge}{\Float}{
 *	        Rate at which the implementation tries to replace the current path
 *	        with a completely new one. Usually, there is little need to change
 *	        this. \default{0.3}
 *	    }
 *	    \parameter{lightImage}{\Boolean}{Include sampling strategies that connect
 *	        paths traced from emitters directly to the camera? (i.e. what \pluginref{ptracer} does)
 *  	    This improves the effectiveness of bidirectional path tracing
 *	        but severely increases the local and remote communication
 *	        overhead, since large \emph{light images} must be transferred between threads
 *	        or over the network. See the text below for a more detailed explanation.
 *	        \default{include these strategies, i.e. \code{true}}
 *      }
 *      \parameter{acceptanceMap}{\Boolean}{
 *          Output acceptance map instead of the rendered image. 
 *          \default{{\footnotesize\code{false}}}
 *      }
 *      \parameter{fixEmitterPath}{\Boolean}{
 *          When using \code{mmlt} setting this flag to true will
 *          only mutate the sensor subpath during the second stage, unless pure emitter. Only 
 *          compatible with \code{mmlt}. \default{{\footnotesize\code{false}}}
 *      }
 *      \parameter{timidAfterLarge}{\Boolean}{
 *          To perform or not, a second stage after a large step. Not compatible
 *          with \code{mira}. \default{{\footnotesize\code{false}}}
 *      }
 *      \parameter{useMixture}{\Boolean}{
 *          Use MH with an equal weight mixture of both stage instead of DR. 
 *          \default{{\footnotesize\code{false}}}
 *      }
 *	    \parameter{sigma}{\Float}{
 *	        If using a Gaussian transition kernel, this correspond to its base 
 *          standard deviation\default{1./64.}
 *	    }
 *	    \parameter{scaleSecond}{\Float}{
 *	        If using a Gaussian transition kernel during the second, this correspond 
 *          the factor applied to sigma to scale it down \default{0.1}
 *	    }
 * }
 * Delayed Rejection Metropolis Light Transport (PSSMLT) is a rendering
 * technique developed by Rioux et al. which is
 * based on Delayed Rejection Monte Carlo (MCMC) integration and Primary 
 * Sample Metropolis Light Transport.
 * 
 * In contrast to the usual methods based on the Metropolis-Hasting algorithm
 * DRMLT uses a two-stage hierarchy of proposal in a single iteration. The goal 
 * of this approach is to fall back on better suited transition kernel in an 
 * informed manner. The resulting is a method that is abe to balance between global
 * exploration and local exploitation.
 * 
 * The DRMLT implementation in Mitsuba implement the \emph{bold-then-timid} approach 
 * discussed in the paper, the \emph{cheap-then-expensive} one being implemented
 * in Li & al (20015) differentiable path tracer DPT. It can operate on top of
 * either a simple unidirectional volumetric path tracer, a fully-fledged bidirectional path
 * tracer with  multiple importance sampling or a multiplexed path tracer, and this 
 * choice is controlled by the \code{technique} flag.
 * 
 * Also, this implementation support Tierney & Mira (1999) original DR , Green & Mira 
 * (2001) generalized one and the pairwise orbital approach to the original framefork.
 * This choice is controlled by the \code{type} flag.
 * 
 * To control the scaling of the second stage mutation size with respect to the first stage,
 * you can use the \code{scaleSecond} parameter.
 * 
 * To perform a second stage after a large step mutation, you can use the \code{timidAfterLarge} 
 * flag.However, doing so is not compatible with mira's approach.
 * 
 * Instead of outputing the rendered image, you can choose to output the acceptance map 
 * of the run using the \code{acceptanceMap} flag. the R channel correspond to the first 
 * stage and G to the second stage accepted mutations.
 * 
 * When using the multiplexed application of DRMLT there's a few caveats to take into 
 * consideration. First, the path sampling technique is kept constant between each 
 * large step along a chain. In our experimentation, doing so does not affect the
 * efficiency of the mmlt while making it easier to do two stage mutations.
 * Secondly, only the sensor subpath is mutated during the second stage, unless the 
 * path is a pure emitter one where we mutate the full path. The control of this feature
 * is done with the \code{fixEmitterPath} flag. Setting this flag to \code{true} with any other
 * technique will result in failure to run.
 * 
 * Finallly, you can compare DRMLT against an equal weight mixture of both stage using 
 * MH by setting the \code{useMixture} flag to \code{True}.
 */


class DRMLT : public Integrator {
public:
    DRMLT(const Properties &props) : Integrator(props) {
        /* Note: a bunch of the parameters below are not publicly exposed,
           because there is really little sense for most users to ever change them. */
        
        
        /* ==================================================================== */
        /*                           General                                    */
        /* ==================================================================== */

        /* If set to <tt>bdpt</tt>, the MLT algorithm runs on top of a
            bidirectional path tracer with multiple importance sampling.
            If set to <tt>mmlt</tt>, the MLT algorithm runs on top of a
            multiplexed path sampler
            If set to <tt>path</tt>, the MLT algorithm runs on top of a
            to a basic path tracer. */
        m_config.technique = [&]() -> PathSampler::ETechnique {
            auto integrator = props.getString("technique");
            if (integrator == "path") {
                return PathSampler::EUnidirectional;
            } else if (integrator == "bdpt") {
                return PathSampler::EBidirectional;
            } else if (integrator == "mmlt") {
                return PathSampler::EMMLT;
            } else {
                Log(EError, "Unknown technique type");
            }
        }();

        
        /* Longest visualized path length (<tt>-1</tt>=infinite).
            A value of <tt>1</tt> will visualize only directly visible light
            sources. <tt>2</tt> will lead to single-bounce (direct-only)
            illumination, and so on. */
        m_config.maxDepth = props.getInteger("maxDepth", -1);

        if (m_config.technique == PathSampler::EMMLT && m_config.maxDepth == -1) {
            Log(EError, "Impossible to use MMLT with no max depth");
        }

        /* Depth to begin using russian roulette (set to -1 to disable) */
        m_config.rrDepth = props.getInteger("rrDepth", 5);


        /* Should an optimized direct illumination sampling strategy be used
            for s=1 paths? (as opposed to plain emission sampling). Usually
            a good idea. Note that this setting only applies when the
            bidirectional path tracer is used internally. The optimization
            affects all paths, not just the ones contributing direct illumination,
            hence it is completely unrelated to the <tt>separateDirect</tt>
            parameter. */
        m_config.directSampling = props.getBoolean("directSampling", true);
        if (m_config.technique == PathSampler::EMMLT) {
            m_config.directSampling = false; 
        }

        /* This parameter can be used to specify the samples per pixel used to
            render the direct component. Should be a power of two (otherwise, it will
            be rounded to the next one). When set to zero or less, the
            direct illumination component will be hidden, which is useful
            for analyzing the component rendered by MLT. When set to -1,
            DRMLT will handle direct illumination as well */
        m_config.directSamples = props.getInteger("directSamples", 16);

        /*Separate direct illumination computation*/
        m_config.separateDirect = m_config.directSamples >= 0;

        /* Number of samples used to estimate the total luminance
           received by the sensor's sensor */
        m_config.luminanceSamples = props.getInteger("luminanceSamples", 100000);
        
        /* Probability of creating large mutations in the [Kelemen et. al]
           MLT variant. The default is 0.3. */
        m_config.pLarge = props.getFloat("pLarge", 0.3f);
    
        /* Specifies the number of parallel work units required for
            multithreaded and network rendering. When set to <tt>-1</tt>, the
            amount will default to four times the number of cores. Note that
            every additional work unit entails a significant amount of
            communication overhead (a full-sized floating put image must be
            transmitted), hence it is important to set this value as low as
            possible, while ensuring that there are enough units to keep all
            workers busy. */
        m_config.workUnits = props.getInteger("workUnits", -1);

        /* Should the multiple importance sampling-based weight computation by
            Kelemen et al. be used? Otherwise, the implementation falls back
            to the 'use of expectations' technique from Veach-style MLT. */
        m_config.kelemenStyleWeights = props.getBoolean("kelemenStyleWeights", true);
        if (m_config.technique == PathSampler::EMMLT) {
            m_config.kelemenStyleWeights = false;
        }

        /* This setting can be very useful to reduce noise in dark regions
            of the image: it activates two-stage MLT, where a nested MLT renderer
            first creates a tiny version of the output image. In a second pass,
            the full version is then rendered, while making use of information
            about the image-space luminance distribution found in the first
            pass. Two-stage MLT is very useful in making the noise characteristics
            more uniform over time image -- specifically, since MLT tends to get
            stuck in very bright regions at the cost of the remainder of the image.*/
        m_config.twoStage = props.getBoolean("twoStage", false);

        

        /* Used internally to let the nested rendering process of a
            two-stage MLT approach know that it is running the first stage */
        m_config.firstStage = props.getBoolean("firstStage", false);

        /* When running two-stage MLT, this parameter determines the size
            of the downsampled image created in the first pass (i.e. setting this
            to 16 means that the horizontal/vertical resolution will be 16 times
            lower). Usually, it's fine to leave this parameter unchanged. When
            the two-stage process introduces noisy halos around very bright image
            regions, it can be set to a lower value */
        m_config.firstStageSizeReduction = props.getInteger(
            "firstStageSizeReduction", 16);

        /* Stop MLT after X seconds -- useful for equal-time comparisons */
        m_config.timeout = props.getInteger("timeout", 0);
        
        /* Use precomputed normalization factor */
        m_config.averageLuminance = props.getFloat("averageLuminance", -1.0f);

        m_config.lightImage = props.getBoolean("lightImage", true);


        /* ==================================================================== */
        /*                           DRMLT Specifics                            */
        /* ==================================================================== */

        /* Select between different delayed rejection implementation,
            <tt>green</tt> for Green & Mira (2001)
            <tt>mira</tt> for Tierney & Mira (1999)
            <tt>orbital</tt> for Orbital */
        m_config.type = [&]() -> DRMLTConfiguration::EType {
            auto implementation = props.getString("type");
            if(implementation == "green") {
                return DRMLTConfiguration::EGreen;
            } else if (implementation == "mira") {
                return DRMLTConfiguration::EMira;
            } else if (implementation == "mirasym" || implementation == "orbital") {
                return DRMLTConfiguration::EOrbital;
            } else {
                Log(EError, "Unknown implementation type");
            }
        }();

        /* if <tt>true</tt> output the acceptance heatmap for DR stages*/
        m_config.acceptanceMap = props.getBoolean("acceptanceMap", false);
        
        /* if <tt>true</tt> perform a second stage after a large step */
        m_config.timidAfterLarge = props.getBoolean("timidAfterLarge", false);
        
        /* if <tt>true</tt> fix the emmiter subpath at the second stage*/
        m_config.fixEmitterPath = props.getBoolean("fixEmitterPath", false);
        if (m_config.technique != PathSampler::EMMLT) {
            if (m_config.fixEmitterPath) {
                SLog(EError, "Impossible to use fixEmitterPath without MMLT");
            }
        }

        /* if <tt>true</tt>  use a mixture instead of two stages*/
        m_config.useMixture = props.getBoolean("useMixture", false);
        
        /* Gaussian Kernel standard deviation */
        m_config.sigma = props.getFloat("sigma", 1.0f / 64.0f);
        
        /* scale factor of the second stage kernel*/
        m_config.scaleSecond = props.getFloat("scaleSecond", 0.1);
        if (m_config.scaleSecond > 1.0) {
            SLog(EError, "scaleSecond is bigger than the first stage");
        }

    }

    /** 
     * Unserialize from a binary data stream
    */
    DRMLT(Stream *stream, InstanceManager *manager)
        : Integrator(stream, manager) {
        m_config = DRMLTConfiguration(stream);
        configure();
    }

    virtual ~DRMLT() {}

    void serialize(Stream *stream, InstanceManager *manager) const {
        Integrator::serialize(stream, manager);
        m_config.serialize(stream);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
                    const RenderJob *job, int sceneResID, int sensorResID,
                    int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                                sensorResID, samplerResID);

        ref<const Sensor> sensor = scene->getSensor();

        if (scene->getSubsurfaceIntegrators().size() > 0)
            Log(EError, "Subsurface integrators are not supported by MLT!");

        if (sensor->getSampler()->getClass()->getName() != "IndependentSampler")
            Log(EError, "Metropolis light transport requires the independent sampler");

        return true;
    }

    void cancel() {
        ref<RenderJob> nested = m_nestedJob;
        if (nested)
            nested->cancel();
        Scheduler::getInstance()->cancel(m_process);
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
                int sceneResID, int sensorResID, int samplerResID) {
        ref<Scheduler> scheduler = Scheduler::getInstance();
        ref<Sensor> sensor = scene->getSensor();
        ref<Sampler> sampler = sensor->getSampler();
        const Film *film = sensor->getFilm();
        size_t nCores = scheduler->getCoreCount();
        size_t sampleCount = sampler->getSampleCount();
        m_config.importanceMap = NULL;

        /* Avoid dumping black images in bootstrap stage */
        scene->setImageDump(false);

        if (m_config.twoStage && !m_config.firstStage) {
            Log(EInfo, "Executing first MLT stage");
            ref<Timer> timer = new Timer();
            Assert(m_config.firstStageSizeReduction > 0);
            m_config.importanceMap = BidirectionalUtils::mltLuminancePass(
                scene, sceneResID, queue, m_config.firstStageSizeReduction,
                m_nestedJob);
            if (!m_config.importanceMap) {
                Log(EWarn, "First-stage MLT process failed!");
                return false;
            }
            Log(EInfo, "First MLT stage took %i ms", timer->getMilliseconds());
        }

        bool nested = m_config.twoStage && m_config.firstStage;

        Vector2i cropSize = film->getCropSize();
        Assert(cropSize.x > 0 && cropSize.y > 0);
        Log(EInfo, "Starting %srender job (%ix%i, "
            SIZE_T_FMT
            " %s, "
            SSE_STR
            ", approx. "
            SIZE_T_FMT
            " mutations/pixel) ..",
            nested ? "nested " : "", cropSize.x, cropSize.y,
            nCores, nCores == 1 ? "core" : "cores", sampleCount);

        auto desiredMutationsPerWorkUnit = [&]() -> size_t {
            if (m_config.technique == PathSampler::EUnidirectional) {
                return 200000;
            } else if(m_config.technique == PathSampler::EBidirectional) {
                return 100000;
            } else if(m_config.technique == PathSampler::EMMLT) {
                return 100000;
            } else {
                SLog(EError, "Unknown path sampler");
            }
        }();

        if (m_config.workUnits <= 0) {
            const size_t cropArea  = (size_t) cropSize.x * cropSize.y;
            const size_t workUnits = ((desiredMutationsPerWorkUnit - 1) +
                (cropArea * sampleCount)) / desiredMutationsPerWorkUnit;
            Assert(workUnits <= (size_t) std::numeric_limits<int>::max());
            m_config.workUnits = (int) std::max(workUnits, (size_t) 1);
        }

        size_t luminanceSamples = m_config.luminanceSamples;
        {
            int times = 10;
            if (m_config.technique == PathSampler::EMMLT) {
                times = 50; // Very conservative initialization
            }
            if (luminanceSamples < (size_t) m_config.workUnits * times) {
                luminanceSamples = (size_t) m_config.workUnits * times;
                Log(EWarn, "Warning: increasing number of luminance samples to "
                    SIZE_T_FMT,
                    luminanceSamples);
            }
        }

        if (m_config.technique == PathSampler::EMMLT) {
            /* In this case we will increase the number of samples
            per max depth (to be consistent to PBRT implementation) */
            SAssert(m_config.maxDepth > 0); // TODO: Depth = 0 is impossible, right?
            luminanceSamples *= (m_config.maxDepth); // Ok
        }

        m_config.nMutations = (cropSize.x * cropSize.y *
            sampleCount) / m_config.workUnits;

        ref<Bitmap> directImage;
        if (m_config.separateDirect && m_config.directSamples > 0 && !nested) {
            directImage = BidirectionalUtils::renderDirectComponent(scene,
                                                                    sceneResID,
                                                                    sensorResID,
                                                                    queue,
                                                                    job,
                                                                    m_config.directSamples);
            if (directImage == NULL)
                return false;
        }

        /* Multi-threaded initialization */
        auto rplSamplers = [&]() -> ref_vector<ReplayableSampler> {
            ref_vector<ReplayableSampler> tmp;
            ref<Random> random = new Random();
            for(int i = 0; i < nCores; i++) tmp.push_back(new ReplayableSampler(random));
            return tmp;
        }();

        std::vector<PathSeed> pathSeeds;
        {
            int nCoresInitialization = nCores;
            if(m_config.technique != PathSampler::EMMLT) {
                nCoresInitialization = 1;
            }

            /* Reducing the number of samples/work united requested as this work
            will be done multiple time by each available thread */
            luminanceSamples = luminanceSamples / nCoresInitialization;
            m_config.workUnits = m_config.workUnits / nCoresInitialization;
            m_config.luminance = 0.0;

            std::mutex merging_mutex;
            BlockScheduler initializeChains(nCoresInitialization, nCoresInitialization, 1);
            initializeChains.run([&](int tileID, int threadID) {
                ref<PathSampler> pathSampler = new PathSampler(m_config.technique,
                                                                scene,
                                                                rplSamplers[threadID].get(),
                                                                rplSamplers[threadID].get(),
                                                                rplSamplers[threadID].get(),
                                                                m_config.maxDepth,
                                                                m_config.rrDepth,
                                                                m_config.separateDirect,
                                                                m_config.directSampling,
                                                                m_config.lightImage);
                std::vector<PathSeed> pathSeedsLocal;
                Float luminance = pathSampler->generateSeeds(luminanceSamples,
                                                                m_config.workUnits, 
                                                                false, 
                                                                m_config.importanceMap,
                                                                pathSeedsLocal);
                /* Merging the information to the main thread */
                {
                    std::lock_guard<std::mutex> lk(merging_mutex);
                    m_config.luminance += luminance;
                    for(auto e: pathSeedsLocal) {
                        e.sampler_id = threadID;
                        pathSeeds.push_back(e);
                    }
                }
            });

            /* Update the number of work unit, luminance */
            luminanceSamples = luminanceSamples * nCoresInitialization;
            m_config.workUnits = m_config.workUnits * nCoresInitialization;
            m_config.luminance /= nCoresInitialization; // Monte carlo estimator
            SLog(EInfo, "Normalization factor computed: %lf", m_config.luminance);
        }


        /* First + second stages heatmap */
        if (m_config.acceptanceMap) {
            m_config.luminance = 1.0f;
        }

        /* Precomputed normalization factor for better comparisons */
        else if (m_config.averageLuminance != -1.0f) {
            m_config.luminance = m_config.averageLuminance;
            Log(EWarn, "Replacing average luminance with the precomputed one = %f", m_config.luminance);
        }

        if (!nested)
            m_config.dump();

        /* Now ready to dump partial images */
        scene->setImageDump(true);

		/* Create the process */
        ref<DRMLTProcess> process = new DRMLTProcess(job, 
                                                        queue,
                                                        m_config, 
                                                        directImage, 
                                                        pathSeeds, 
                                                        rplSamplers);

        /* Create a sampler instance for each worker */
        ref<DRMLTSampler> mltSampler;
        switch (m_config.type) {
            case DRMLTConfiguration::EGreen:
                mltSampler = new GreenDRMLTSampler(m_config);
                break;
            case DRMLTConfiguration::EMira:
                mltSampler = new MiraDRMLTSampler(m_config);
                break;
            case DRMLTConfiguration::EOrbital:
                mltSampler = new OrbitalDRMLTSampler(m_config);
                break;
            default:
                SLog(EError, "Invalid sampler type");
                break;
        }
        std::vector<SerializableObject *> mltSamplers(scheduler->getCoreCount());
        for (size_t i = 0; i < mltSamplers.size(); ++i) {
            ref<Sampler> clonedSampler = mltSampler->clone();
            clonedSampler->incRef();
            mltSamplers[i] = clonedSampler.get();
        }
        int mltSamplerResID = scheduler->registerMultiResource(mltSamplers);
        for (size_t i = 0; i < scheduler->getCoreCount(); ++i)
            mltSamplers[i]->decRef();

        process->bindResource("scene", sceneResID);
        process->bindResource("sensor", sensorResID);
        process->bindResource("sampler", mltSamplerResID);

        m_process = process;
        scheduler->schedule(process);
        scheduler->wait(process);
        m_process = NULL;
        process->develop();

        return process->getReturnStatus() == ParallelProcess::ESuccess;
    }

    MTS_DECLARE_CLASS()
private:
    ref<ParallelProcess> m_process;
    ref<RenderJob> m_nestedJob;
    DRMLTConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(DRMLT, false, Integrator)
MTS_EXPORT_PLUGIN(DRMLT, "Delayed Rejection MLT");
MTS_NAMESPACE_END
