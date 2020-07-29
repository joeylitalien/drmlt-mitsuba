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
#include <mitsuba/bidir/path.h>
#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/bidir/rsampler.h>
#include "drmlt_proc.h"
#include "drmlt_sampler.h"

#include "../pssmlt_utils.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Rendering Statistics                         */
/* ==================================================================== */

StatsCounter firstLevelRatio("Delayed Rejection MLT", 
    "Accepted 1st-stage mutations", EPercentage);
StatsCounter largeStepRatio("Delayed Rejection MLT", 
    "Accepted large mutations in the 1st-stage mutations", EPercentage);
StatsCounter boldStepRatio("Delayed Rejection MLT", 
    "Accepted bold mutation in the 1st-stage mutations", EPercentage);
StatsCounter secondLevelRatio("Delayed Rejection MLT", 
    "Accepted 2nd-stage mutations", EPercentage);
StatsCounter secondLevelLargeRatio("Delayed Rejection MLT", 
    "Accepted 2nd-stage mutations after large mutation", EPercentage);
StatsCounter secondLevelBoldRatio("Delayed Rejection MLT", 
    "Accepted 2nd-stage mutations after bold mutation", EPercentage);
StatsCounter acceptanceRate("Delayed Rejection MLT", 
    "Overall acceptance rate", EPercentage);
StatsCounter forcedAcceptance("Delayed Rejection MLT", 
    "Number of forced acceptances");

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

class DRMLTRenderer : public WorkProcessor {
public:
    DRMLTRenderer(const DRMLTConfiguration &conf, const ref_vector<ReplayableSampler> rplSamplers)
        : m_config(conf), m_rplSamplers(rplSamplers) {
    }

    DRMLTRenderer(Stream *stream, InstanceManager *manager)
        : WorkProcessor(stream, manager) {
        m_config = DRMLTConfiguration(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        m_config.serialize(stream);
    }

    ref<WorkUnit> createWorkUnit() const {
        return new SeedWorkUnit();
    }

    ref<WorkResult> createWorkResult() const {
        /* Hack to detect box filter */
        Float radius = m_film->getReconstructionFilter()->getRadius();
        if (m_config.acceptanceMap && (radius - 0.500010f > 1e-6)) {
            Log(EError, "Box filter required for acceptance map!");
        }
        return new ImageBlock(Bitmap::ESpectrum,
                                m_film->getCropSize(), m_film->getReconstructionFilter());
    }

    void prepare() {
        Scene *scene = static_cast<Scene *>(getResource("scene"));
        switch (m_config.type) {
            /* Green and Mira (2001) */
            case DRMLTConfiguration::EGreen:
                m_origSampler = static_cast<GreenDRMLTSampler *>(getResource("sampler"));
                m_sensorSampler = new GreenDRMLTSampler(m_origSampler);
                m_emitterSampler = new GreenDRMLTSampler(m_origSampler);
                m_directSampler = new GreenDRMLTSampler(m_origSampler);
                break;
            /* Tierney and Mira (1999) */
            case DRMLTConfiguration::EMira:
                m_origSampler = static_cast<MiraDRMLTSampler *>(getResource("sampler"));
                m_sensorSampler = new MiraDRMLTSampler(m_origSampler);
                m_emitterSampler = new MiraDRMLTSampler(m_origSampler);
                m_directSampler = new MiraDRMLTSampler(m_origSampler);
                break;
            /* Pairwise Orbital */
            case DRMLTConfiguration::EOrbital:
                m_origSampler = static_cast<OrbitalDRMLTSampler *>(getResource("sampler"));
                m_sensorSampler = new OrbitalDRMLTSampler(m_origSampler);
                m_emitterSampler = new OrbitalDRMLTSampler(m_origSampler);
                m_directSampler = new OrbitalDRMLTSampler(m_origSampler);
                break;
            default:
                SLog(EError, "Invalid sampler type");
                break;
        }
        m_sensor = static_cast<Sensor *>(getResource("sensor"));
        m_scene = new Scene(scene);
        m_film = m_sensor->getFilm();
        m_scene->setSensor(m_sensor);
        m_scene->setSampler(m_origSampler);
        m_scene->removeSensor(scene->getSensor());
        m_scene->addSensor(m_sensor);
        m_scene->setSensor(m_sensor);
        m_scene->wakeup(NULL, m_resources);
        m_scene->initializeBidirectional();

		/* Clone all the rplSampler to avoid any race condition between the threads */
        ref_vector<ReplayableSampler> cloned_rplSamplers;
        for(auto i = 0; i < m_rplSamplers.size(); i++) {
            auto tmp = m_rplSamplers[i]->clone();
            auto cast = (ReplayableSampler*)tmp.get();
            cloned_rplSamplers.push_back(cast);
        }
        m_rplSamplers = cloned_rplSamplers;

        /* when using EMMLT, use the direct sampler for the strategy component */
        if (m_config.technique == PathSampler::EMMLT) {
            /* Strategy are kept fixed -> set direct samples stages to identity */
            m_directSampler->setStagesToIdentity();
            if (m_config.fixEmitterPath) {
                /* if fixed emitter subpath, set second stage as identity (unless 
                full emitter path) */
                m_emitterSampler->handleLightTracing();
            }
        }

        m_pathSampler = new PathSampler(m_config.technique, 
                                        m_scene,
                                        m_emitterSampler, 
                                        m_sensorSampler, 
                                        m_directSampler, 
                                        m_config.maxDepth,
                                        m_config.rrDepth, 
                                        m_config.separateDirect, 
                                        m_config.directSampling, 
                                        m_config.lightImage);

    }

    /** 
     * Called in process if useMixture is true. Compute th results for the mixture case 
     * in the EMMLT application. Here the mixture is made of the same transition kernel 
     * as the 2 stages of delayed rejection with a mixture weight of 0.5.
    */
    void processMixture(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {

        ImageBlock *result = static_cast<ImageBlock *>(workResult);
        const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
        const PathSeed &seed = wu->getSeed();
        auto current = std::make_unique<SplatList>();
        auto proposed = std::make_unique<SplatList>();

        /* Useful lambdas */

        /* Good old coin toss */
        auto flipCoin = [](Float x, Random *random) {
            SAssert(x >= 0 && x <= 1.0);
            return (x >= 1) || (random->nextFloat() < x);
        };
        
        /* Metropolis-Hasting acceptance formula */
        auto metropolisClamp = [](Float x) { return std::min((Float) 1.0f, x); };
        
        /* Check if luminance is invalid */
        auto isInvalid = [](Float x) { return std::isnan(x) || std::isinf(x)  || x < 0; };
        
        /* Splat a SplatList with appropriate weight */
        auto splat = [&](SplatList *splat, Float weight) {
            if (weight > 0) {
                for (size_t k = 0; k < splat->size(); ++k) {
                    Spectrum value = splat->getValue(k) * weight;
                    if (value.isValid())
                        result->put(splat->getPosition(k), &value[0]);
                    else
                        SLog(EWarn, "Invalid splat");
                }
            }
        };

        /* Is -1 by default, a certain value for MMLT */
        int depth = seed.depth;
        {
            /* For simplicity, we replace the usual lazy evaluation
            procedure by a full state mutation. All state are assumed to
            be of fixed dimension. We set the dimension of each sampler
            as the worst case they can be given the technique, scene 
            parameters and maximum depth. */
            auto maxDim = findMaxDimensions(m_scene, m_config.maxDepth, m_config.rrDepth, depth,
                                            m_config.technique, m_config.directSampling);
            m_sensorSampler->setMaxDim(maxDim.sensor);   
            m_emitterSampler->setMaxDim(maxDim.emitter); 
            m_directSampler->setMaxDim(maxDim.direct); 
        }

        m_emitterSampler->reset();
        m_sensorSampler->reset();
        m_directSampler->reset();
        m_sensorSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());
        m_emitterSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());
        m_directSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());

        /* Generate the initial sample by replaying the seeding random
        number stream at the appropriate position. Afterwards, revert
        back to this worker's own source of random numbers */
        m_rplSamplers[seed.sampler_id]->setSampleIndex(seed.sampleIndex);
        m_sensorSampler->setReplay(true);
        m_emitterSampler->setReplay(true);
        m_directSampler->setReplay(true);
        m_pathSampler->sampleSplats(Point2i(-1), *current, depth);
        m_sensorSampler->setReplay(false);
        m_emitterSampler->setReplay(false);
        m_directSampler->setReplay(false);
        result->clear();

        ref<Random> random = m_origSampler->getRandom();
        m_sensorSampler->setRandom(random);
        m_emitterSampler->setRandom(random);
        m_directSampler->setRandom(random);
        m_rplSamplers[seed.sampler_id]->updateSampleIndex(m_rplSamplers[seed.sampler_id]->getSampleIndex()
                                        + m_sensorSampler->getSampleIndex()
                                        + m_emitterSampler->getSampleIndex()
                                        + m_directSampler->getSampleIndex());
        /* set replayed state as current */
        m_sensorSampler->accept(true);
        m_emitterSampler->accept(true);
        m_directSampler->accept(true);

        /* ensure that current is filled up to the maximal dimensionality */
        m_sensorSampler->fillReplay();
        m_emitterSampler->fillReplay();
        m_directSampler->fillReplay();

        /* Sanity check -- the luminance should match the one from
        the warmup phase - an error here would indicate inconsistencies
        regarding the use of random numbers during sample generation */
        if (std::abs((current->luminance - seed.luminance)
                        / seed.luminance) > Epsilon)
            Log(EError, "Error when reconstructing a seed path: luminance "
                        "= %f, but expected luminance = %f", current->luminance, seed.luminance);

        current->normalize(m_config.importanceMap);
        /// MLT main loop
        ref<Timer> timer = new Timer();
        for (uint64_t mutationCtr = 0; mutationCtr < m_config.nMutations && !stop; ) {
            if (wu->getTimeout() > 0 && (mutationCtr % 8192) == 0
                && (int) timer->getMilliseconds() > wu->getTimeout())
                break;

            if (isInvalid(current->luminance)) {
                Log(EWarn, "Current invalid");
            }

            /* Reset acceptance variables of previous iteration  */
            Float a = 0.f;
            bool accept = false;

            /* Check if large step, if so set it */
            bool largeStep = random->nextFloat() < m_config.pLarge;
            m_sensorSampler->setLargeStep(largeStep);
            m_emitterSampler->setLargeStep(largeStep);
            m_directSampler->setLargeStep(largeStep);

            /* Generate/Trace first stage proposal */
            m_pathSampler->sampleSplats(Point2i(-1), *proposed, depth);
            proposed->normalize(m_config.importanceMap);
            ++mutationCtr;

            /* Accept through regular MH  */
            if (isInvalid(proposed->luminance)) {
                Log(EWarn, "Encountered a sample with luminance = %f, ignoring! at first stage",
                    proposed->luminance);
                a = 0;
                accept = false;
            } else {
                a = metropolisClamp(proposed->luminance / current->luminance);
                accept = flipCoin(a, random);
            }

            /* Mixture case: proceed to second stage with prob 0.5 and if not large step */
            bool doSecond = false;
            if (!largeStep) {
                doSecond = flipCoin(0.5, random);
            }

            /* Bump up proposal to second stage */
            if (doSecond) {
                /* Set samplers transition kernel to second stage */
                m_sensorSampler->nextStage();
                m_directSampler->nextStage();

                /* EMMLT: Light tracing scenario, switch kernel depending if pure emitter path */
                if(m_config.fixEmitterPath) {
                    m_emitterSampler->nextStage(current->t == 1);
                } else {
                    m_emitterSampler->nextStage();
                }
                /* Generate/Trace second stage proposal (rewrite proposed) */
                m_pathSampler->sampleSplats(Point2i(-1), *proposed, depth);
                proposed->normalize(m_config.importanceMap);

                /* Mixture case: Accept through regular MH (rewrite a and accept) */
                if (isInvalid(proposed->luminance)) {
                    a = 0;
                    accept = false;
                } else {
                    a = metropolisClamp(proposed->luminance / current->luminance);
                    accept = flipCoin(a, random);
                }
            }

            /* Splat states according to their relative weight.*/
            {
                Float proposedWeight = a;
                Float currentWeight = 1.0f - proposedWeight;
                splat(current.get(), currentWeight);
                splat(proposed.get(), proposedWeight);
            }

            /* Either accept proposed state */
            if (accept) {
                proposed.swap(current);
                m_sensorSampler->accept(!doSecond);
                m_emitterSampler->accept(!doSecond);
                m_directSampler->accept(!doSecond);

                /* Stats computation */
                acceptanceRate.incrementBase(1);
                ++acceptanceRate;
                if (!doSecond) {
                    firstLevelRatio.incrementBase(1);
                    ++firstLevelRatio;
                    if (largeStep) {
                        largeStepRatio.incrementBase(1);
                        ++largeStepRatio;
                    } else {
                        boldStepRatio.incrementBase(1);
                        ++boldStepRatio;
                    }
                } else {
                    secondLevelRatio.incrementBase(1);
                    ++secondLevelRatio;
                }
            } 
            /* Or reject proposed state */
            else {
                m_sensorSampler->reject();
                m_emitterSampler->reject();
                m_directSampler->reject();

                /* Stats computation */
                acceptanceRate.incrementBase(1);
                if(!doSecond) {
                    firstLevelRatio.incrementBase(1);
                    if (largeStep) {
                        largeStepRatio.incrementBase(1);
                    } else {
                        boldStepRatio.incrementBase(1);
                    }
                } else {
                    secondLevelRatio.incrementBase(1);
                }
            }
        }
    }

    /** 
     * Core of the Delayed Rejection Metropolis Light Transport DRMLT (2020).
     * in case of useMixture being true, it will call the processMixture function instead
    */
    void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
        
        /* If using regular MH with the mixture composed of both DR stages */
        if(m_config.useMixture) {
            processMixture(workUnit, workResult, stop);
            return;
        }

        ImageBlock *result = static_cast<ImageBlock *>(workResult);
        const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
        const PathSeed &seed = wu->getSeed();
        auto current = std::make_unique<SplatList>();

        /* First and second stage acceptance colors (for visualization)  */
        bool drawAcceptanceMap = m_config.acceptanceMap;
        const bool VERBOSE = false;
        std::pair<Color3, Color3> stageColors = {Color3(1,0,0), Color3(0,1,0)};

        /* Vector of proposed splat lists for delayed rejection  */
        std::pair<std::unique_ptr<SplatList>, std::unique_ptr<SplatList>> proposed;
        proposed = {std::make_unique<SplatList>(), std::make_unique<SplatList>()};

        /* Splatlist for reverse move y* */
        auto reverse = std::make_unique<SplatList>();

        /* Proposal acceptances */
        const bool splatReverse = m_config.type != DRMLTConfiguration::EMira && false;
        const Float secondProb = 1.0;


        /* Useful lambdas */

        /* Good old coin toss */
        auto flipCoin = [](Float x, Random *random) {
            SAssert(x >= 0 && x <= 1.0);
            return (x >= 1) || (random->nextFloat() < x);
        };

        /* Metropolis-Hasting acceptance formula */
        auto metropolisClamp = [](Float x) { return std::min((Float) 1.0f, x); };
        
        /* Check if luminance is invalid */
        auto isInvalid = [](Float x) { return std::isnan(x) || std::isinf(x)  || x <= 0; };

        /* Splat a SplatList with appropriate weight */
        auto splat = [&](SplatList *splat, Float weight) {
            if (!drawAcceptanceMap && weight > 0) {
                for (size_t k = 0; k < splat->size(); ++k) {
                    Spectrum value = splat->getValue(k) * weight;
                    if (value.isValid())
                        result->put(splat->getPosition(k), &value[0]);
                    else
                        SLog(EWarn, "Invalid splat");
                }
            }
        };

        auto splatAcceptanceOnly = [&](SplatList *splat, int stage) {
            if (drawAcceptanceMap) {
                for (size_t k = 0; k < splat->size(); ++k) {
                    Spectrum value = stage == 0 ? stageColors.first : stageColors.second;
                    result->put(splat->getPosition(k), &value[0]);
                }
            }
        };

        /* Is -1 by default, a certain value for MMLT */
        int depth = seed.depth;
        {
            /* For simplicity, we replace the usual lazy evaluation
            procedure by a full state mutation. All state are assumed to
            be of fixed dimension. We set the dimension of each sampler
            as the worst case they can be given the technique, scene 
            parameters and maximum depth. */
            auto maxDim = findMaxDimensions(m_scene, m_config.maxDepth, m_config.rrDepth, depth,
                                            m_config.technique, m_config.directSampling);
            m_sensorSampler->setMaxDim(maxDim.sensor);
            m_emitterSampler->setMaxDim(maxDim.emitter);
            m_directSampler->setMaxDim(maxDim.direct);
        }

        m_emitterSampler->reset();
        m_sensorSampler->reset();
        m_directSampler->reset();
        m_sensorSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());
        m_emitterSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());
        m_directSampler->setRandom(m_rplSamplers[seed.sampler_id]->getRandom());

        /* Generate the initial sample by replaying the seeding random
        number stream at the appropriate position. Afterwards, revert
        back to this worker's own source of random numbers */
        m_rplSamplers[seed.sampler_id]->setSampleIndex(seed.sampleIndex);
        m_sensorSampler->setReplay(true);
        m_emitterSampler->setReplay(true);
        m_directSampler->setReplay(true);
        m_pathSampler->sampleSplats(Point2i(-1), *current, depth);
        m_sensorSampler->setReplay(false);
        m_emitterSampler->setReplay(false);
        m_directSampler->setReplay(false);
        result->clear();

        ref<Random> random = m_origSampler->getRandom();
        m_sensorSampler->setRandom(random);
        m_emitterSampler->setRandom(random);
        m_directSampler->setRandom(random);
        m_rplSamplers[seed.sampler_id]->updateSampleIndex(m_rplSamplers[seed.sampler_id]->getSampleIndex()
                                            + m_sensorSampler->getSampleIndex()
                                            + m_emitterSampler->getSampleIndex()
                                            + m_directSampler->getSampleIndex());

        /* set replayed state as current */
        m_sensorSampler->accept(true);
        m_emitterSampler->accept(true);
        m_directSampler->accept(true);

        /* ensure that current is filled up to the maximal dimensionality */
        m_sensorSampler->fillReplay();
        m_emitterSampler->fillReplay();
        m_directSampler->fillReplay();

        /* Sanity check -- the luminance should match the one from
        the warmup phase - an error here would indicate inconsistencies
        regarding the use of random numbers during sample generation */
        if (std::abs((current->luminance - seed.luminance)
                        / seed.luminance) > Epsilon)
            Log(EError, "Error when reconstructing a seed path: luminance "
                        "= %f, but expected luminance = %f", current->luminance, seed.luminance);

        current->normalize(m_config.importanceMap);

        /* MLT main loop */
        ref<Timer> timer = new Timer();
        for (uint64_t mutationCtr = 0; mutationCtr < m_config.nMutations && !stop; ) {
            if (wu->getTimeout() > 0 && (mutationCtr % 8192) == 0
                && (int) timer->getMilliseconds() > wu->getTimeout())
                break;

            if(isInvalid(current->luminance)) {
                Log(EError, "Current invalid, quit");
            }

            /* Reset acceptance variables of previous iteration  */
            std::pair<Float, Float> a = {0.f, 0.f};
            std::pair<bool, bool> accept = {false, false};
            

            /* Check if large step, if so set it */
            bool largeStep = random->nextFloat() < m_config.pLarge;
            m_sensorSampler->setLargeStep(largeStep);
            m_emitterSampler->setLargeStep(largeStep);
            m_directSampler->setLargeStep(largeStep);

            /* Generate/Trace first stage proposal */
            m_pathSampler->sampleSplats(Point2i(-1), *proposed.first, depth);
            proposed.first->normalize(m_config.importanceMap);
            ++mutationCtr;

            /* Accept through regular MH. See DRMLT (2020), Eq. 5 */
            if (isInvalid(proposed.first->luminance)) {
                a.first = 0;
                accept.first = false;
            } else {
                a.first = metropolisClamp(proposed.first->luminance / current->luminance);
                accept.first = flipCoin(a.first, random);
            }

            /* DR case: Probabilistically proceed to second stage if first is rejected */
            bool doSecond = false;
            doSecond = !accept.first && flipCoin(secondProb, random); // secondProb = 1 so always
            /* Only do second stage after large step if enabled */
            if (!m_config.timidAfterLarge) {
                doSecond = doSecond && !largeStep;
            }

            /* Bump up proposal to second stage */
            if (doSecond) {

                /* Set samplers transition kernel to second stage */
                m_sensorSampler->nextStage();
                m_directSampler->nextStage();

                /* EMMLT: Light tracing scenario, switch kernel depending if pure emitter path */
                if(m_config.fixEmitterPath) {
                    m_emitterSampler->nextStage(current->t == 1);
                } else {
                    m_emitterSampler->nextStage();
                }


                /* Generate/Trace second stage proposal */
                m_pathSampler->sampleSplats(Point2i(-1), *proposed.second, depth);
                proposed.second->normalize(m_config.importanceMap);

                /* DR case: accept with corrected second stage acceptance (depend on type) */
                if (isInvalid(proposed.second->luminance)) {
                    a.second = 0;
                    accept.second = false;
                } else {
                    /**
                     * Green & Mira (2001): use reverse path y^* to remove transition
                     * kernel ratio
                    */
                    if (m_config.type == DRMLTConfiguration::EGreen) {

                        /* Set sampler in reverse mode */
                        m_sensorSampler->setReverse(true);
                        m_emitterSampler->setReverse(true);
                        m_directSampler->setReverse(true);

                        /* Trace reverse path */
                        m_pathSampler->sampleSplats(Point2i(-1), *reverse, depth);
                        reverse->normalize(m_config.importanceMap);

                        /* Compute the fist stage reverse path acceptance ratio. See DRMLT (2020), Eq. 13 */
                        Float aReverse;
                        if (isInvalid(reverse->luminance)) {
                            aReverse = 0.f;
                        } else {
                            aReverse = metropolisClamp(reverse->luminance / proposed.second->luminance);
                        }

                        /* Accept through Green's second stage acceptance formula. See DRMLT (2020), Eq. 14 */
                        if (aReverse == 1) {
                            a.second = 0.f;
                            accept.second = false;
                        } else {
                            auto lumRatio = proposed.second->luminance / current->luminance;
                            a.second = metropolisClamp(lumRatio * (1.f - aReverse) / (1.f - a.first));
                            accept.second = flipCoin(a.second, random);
                        }

                        /* Exit reverse mode */
                        m_sensorSampler->setReverse(false);
                        m_emitterSampler->setReverse(false);
                        m_directSampler->setReverse(false);
                    }
                    /**
                     * Tierney & Mira (1999): Naive approach
                    */
                    else if (m_config.type == DRMLTConfiguration::EMira) {
                        /* compute the first stage reverse acceptance. See DRMLT (2020), Eq. 5 with x = z */
                        auto aReverse = metropolisClamp(proposed.first->luminance / proposed.second->luminance);
                        if (aReverse >= 1) {
                            a.second = 0.0f;
                            accept.second = false;
                        } else {
                            /* Compute transition kernels ratio */
                            Float transitionRatio = largeStep ? 1.0 :
                                                    m_sensorSampler->getTransitionRatio(m_config.pLarge) *
                                                    m_emitterSampler->getTransitionRatio(m_config.pLarge) *
                                                    m_directSampler->getTransitionRatio(m_config.pLarge);
                            SAssert(!largeStep);

                            /* Accept with the naive second stage acceptance formula. See DRMLT (2020), Eq. 7 */
                            if (isInvalid(transitionRatio)) {
                                a.second = 0.0f;
                                accept.second = false;
                            } else {
                                auto lumRatio = proposed.second->luminance / current->luminance;
                                a.second = metropolisClamp(
                                        lumRatio * transitionRatio * (1.0f - aReverse) / (1.0f - a.first));
                                accept.second = flipCoin(a.second, random);
                            }
                        }
                    }
                    /**
                     * Pairwise Orbital: uses Tierney & Mira (1999) with our second stage
                     * pairwise orbital proposal to remove transition kernel ratio.
                    */
                    else if (m_config.type == DRMLTConfiguration::EOrbital) {
                        
                        /* Accept with the simplified second stage acceptance formula. See DRMLT (2020), Eq. 11 */
                        if (proposed.second->luminance < proposed.first->luminance) {
                            a.second = 0.0f;
                            accept.second = false;
                        }
                        else if (proposed.second->luminance >= current->luminance) {
                            a.second = 1.0f;
                            accept.second = true;
                        } else {
                            a.second = (proposed.second->luminance - proposed.first->luminance) /
                                (current->luminance - proposed.first->luminance);
                            accept.second = flipCoin(a.second, random);
                        }
                    } else {
                        SLog(EError, "Invalid algorithm type");
                    }
                }
            }

            /* Compute acceptance weights for splatting. See DRMLT (2020), Fig. 10*/
            Float currentWeight;
            std::pair<Float, Float> proposedWeight;
            proposedWeight.first = a.first;
            proposedWeight.second = (1.0f - a.first) * a.second;
            currentWeight = 1.0f - proposedWeight.first - proposedWeight.second;

            /* Splat current  */
            splat(current.get(), currentWeight);
            /* Splat first stage proposal */
            splat(proposed.first.get(), proposedWeight.first);
            /* Splat second stage proposal */
            splat(proposed.second.get(), proposedWeight.second);

            /* Either accept one of the proposed states */
            if (accept.first || accept.second) {
                bool acceptedFirst = true;
                /* Either first stage was accepted */
                if (accept.first) {
                    proposed.first.swap(current);

                    /* Add contribution to acceptance map */
                    if (!largeStep) {
                        splatAcceptanceOnly(proposed.first.get(), 0);
                    }
                } 
                /* Or irst stage was rejected and the second one accepted */
                else {
                    acceptedFirst = false;
                    proposed.second.swap(current);

                    /* Add contribution to acceptance map */
                    splatAcceptanceOnly(proposed.second.get(), 1);
                }

                m_sensorSampler->accept(accept.first);
                m_emitterSampler->accept(accept.first);
                m_directSampler->accept(accept.first);

                /* Stats computation */
                acceptanceRate.incrementBase(1);
                ++acceptanceRate;
                if (acceptedFirst) {
                    firstLevelRatio.incrementBase(1);
                    ++firstLevelRatio;
                    if (largeStep) {
                        largeStepRatio.incrementBase(1);
                        ++largeStepRatio;
                    } else {
                        boldStepRatio.incrementBase(1);
                        ++boldStepRatio;
                    }
                } else {
                    acceptanceRate.incrementBase(1);
                    firstLevelRatio.incrementBase(1);
                    secondLevelRatio.incrementBase(1);
                    ++secondLevelRatio;
                    if (largeStep) {
                        largeStepRatio.incrementBase(1);
                        secondLevelLargeRatio.incrementBase(1);
                        ++secondLevelLargeRatio;
                    } else {
                        boldStepRatio.incrementBase(1);
                        secondLevelBoldRatio.incrementBase(1);
                        ++secondLevelBoldRatio;
                    }
                }
            } 

            /* Or reject proposed states */
            else {
                m_sensorSampler->reject();
                m_emitterSampler->reject();
                m_directSampler->reject();

                /* Stats computation */
                acceptanceRate.incrementBase(1);
                firstLevelRatio.incrementBase(1);
                if (largeStep) {
                    largeStepRatio.incrementBase(1);
                    if (doSecond) {
                        secondLevelRatio.incrementBase(1);
                        secondLevelLargeRatio.incrementBase(1);
                        acceptanceRate.incrementBase(1);
                    }
                } else {
                    boldStepRatio.incrementBase(1);
                    if (doSecond) {
                        secondLevelRatio.incrementBase(1);
                        secondLevelBoldRatio.incrementBase(1);
                        acceptanceRate.incrementBase(1);
                    }
                }
            }
        }
    }

    ref<WorkProcessor> clone() const {
        return new DRMLTRenderer(m_config, m_rplSamplers);
    }

    MTS_DECLARE_CLASS()
private:
    DRMLTConfiguration m_config;
    ref<Scene> m_scene;
    ref<Sensor> m_sensor;
    ref<Film> m_film;
    ref<PathSampler> m_pathSampler;
    ref<DRMLTSampler> m_origSampler;
    ref<DRMLTSampler> m_sensorSampler;
    ref<DRMLTSampler> m_emitterSampler;
    ref<DRMLTSampler> m_directSampler;
    ref_vector<ReplayableSampler> m_rplSamplers;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

DRMLTProcess::DRMLTProcess(const RenderJob *parent, RenderQueue *queue,
                            const DRMLTConfiguration &conf, const Bitmap *directImage,
                            const std::vector<PathSeed> &seeds, ref_vector<ReplayableSampler> rplSamplers)
                            : m_job(parent), m_queue(queue), m_config(conf),
                            m_progress(NULL), m_seeds(seeds), m_rplSamplers(rplSamplers) {
    m_directImage = directImage;
    m_timeoutTimer = new Timer();
    m_refreshTimer = new Timer();
    m_resultMutex = new Mutex();
    m_resultCounter = 0;
    m_workCounter = 0;
    m_refreshTimeout = 1;
}

ref<WorkProcessor> DRMLTProcess::createWorkProcessor() const {
    return new DRMLTRenderer(m_config, m_rplSamplers);
}

void DRMLTProcess::develop() {
    LockGuard lock(m_resultMutex);
    size_t pixelCount = m_accum->getBitmap()->getPixelCount();
    const Spectrum *accum = (Spectrum *) m_accum->getBitmap()->getData();
    const Spectrum *direct = m_directImage != NULL ?
                             (Spectrum *) m_directImage->getData() : NULL;
    const Float *importanceMap = m_config.importanceMap != NULL ?
                                    m_config.importanceMap->getFloatData() : NULL;
    Spectrum *target = (Spectrum *) m_developBuffer->getData();

    /* Compute the luminance correction factor */
    Float avgLuminance = 0;
    if (importanceMap) {
        for (size_t i = 0; i < pixelCount; ++i)
            avgLuminance += accum[i].getLuminance() * importanceMap[i];
    } else {
        for (size_t i = 0; i < pixelCount; ++i)
            avgLuminance += accum[i].getLuminance();
    }

    avgLuminance /= (Float) pixelCount;
    Float luminanceFactor = m_config.luminance / avgLuminance;

    /* If acceptance map, set factor to 1 (bins) */
    if(m_config.acceptanceMap) {
        luminanceFactor = 1.0;
    }

    for (size_t i = 0; i < pixelCount; ++i) {
        Float correction = luminanceFactor;
        if (importanceMap)
            correction *= importanceMap[i];
        Spectrum value = accum[i] * correction;
        if (direct)
            value += direct[i];
        target[i] = value;
    }
    m_film->setBitmap(m_developBuffer);
    m_refreshTimer->reset();

    m_queue->signalRefresh(m_job);
}

void DRMLTProcess::processResult(const WorkResult *wr, bool cancelled) {
    LockGuard lock(m_resultMutex);
    const ImageBlock *result = static_cast<const ImageBlock *>(wr);
    m_accum->put(result);
    m_progress->update(++m_resultCounter);
    m_refreshTimeout = std::min(2000U, m_refreshTimeout * 2);

    /* Re-develop the entire image every two seconds if partial results are
       visible (e.g. in a graphical user interface). */
    if (m_job->isInteractive() && m_refreshTimer->getMilliseconds() > m_refreshTimeout)
        develop();
}

ParallelProcess::EStatus DRMLTProcess::generateWork(WorkUnit *unit, int worker) {
    int timeout = 0;
    if (m_config.timeout > 0) {
        timeout = static_cast<int>(static_cast<int64_t>(m_config.timeout * 1000) -
            static_cast<int64_t>(m_timeoutTimer->getMilliseconds()));
    }

    if (m_workCounter >= m_config.workUnits || timeout < 0)
        return EFailure;

    SeedWorkUnit *workUnit = static_cast<SeedWorkUnit *>(unit);
    workUnit->setSeed(m_seeds[m_workCounter++]);
    workUnit->setTimeout(timeout);
    return ESuccess;
}

void DRMLTProcess::bindResource(const std::string &name, int id) {
    ParallelProcess::bindResource(name, id);
    if (name == "sensor") {
        m_film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();
        if (m_progress)
            delete m_progress;
        m_progress = new ProgressReporter("Rendering", m_config.workUnits, m_job);
        m_accum = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize());
        m_accum->clear();
        m_developBuffer = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_film->getCropSize());
    }
}

MTS_IMPLEMENT_CLASS_S(DRMLTRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(DRMLTProcess, false, ParallelProcess)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)

MTS_NAMESPACE_END
