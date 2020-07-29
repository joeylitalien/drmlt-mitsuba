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


#if !defined(__DRMLT_SAMPLER_H)
#define __DRMLT_SAMPLER_H

#include <mitsuba/render/sampler.h>
#include <mitsuba/core/random.h>

#include "drmlt.h"
#include "tools/transition.h"

MTS_NAMESPACE_BEGIN

/**
 * Sampler implementation as described in
 * 'Delayed Rejection Metropolis Light Transport 
 * Algorithm' by Rioux et al. This is an abstract 
 * class for all type of DR algoritmn. 
*/
class DRMLTSampler : public Sampler {
public:

    /* ==================================================================== */
    /*                           General                                    */
    /* ==================================================================== */
    
    /// Construct a new MLT sampler
    explicit DRMLTSampler(const DRMLTConfiguration &conf);

    /* Construct a new sampler, which operates on the
	same random number generator as a sampler.*/
    explicit DRMLTSampler(DRMLTSampler *sampler);

    /// Unserialize from a binary data stream
    DRMLTSampler(Stream *stream, InstanceManager *manager);

    /// Setup the internal states
    void configure();

    /// Serialize to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    /// Set whether the current step should be large
    inline void setLargeStep(bool value) { m_largeStep = value; }

    /// Check if the current step is a large step
    inline bool isLargeStep() const { return m_largeStep; }

    /// Retrieve the next component value from the current sample
    virtual Float next1D();

    /// Retrieve the next two component values from the current sample
    virtual Point2 next2D();

    /// Return a string description
    virtual std::string toString() const;

    /// Set replay mode (replay rng from the bootstrap)
    inline void setReplay(bool value) { m_replay = value; }

    /// Set maximum dimension of primary sample vector
    inline void setMaxDim(size_t dim) { m_maxDim = dim; }

    /// Return a primary sample
    virtual Float primarySample(size_t i);

    /// Reset (& start with a uniform mutation)
    void reset();

    /// Accept a mutation
    void accept(bool acceptFirst);

    /// Reject a mutation
    void reject();

    /// Replace the underlying random number generator
    inline void setRandom(Random *random) { m_random = random; }

    /// Return the underlying random number generator
    inline Random *getRandom() { return m_random; }

    /// The following functions do nothing in this implementation
    virtual void advance() {}
    virtual void generate(const Point2i &pos) {}

    /// The following functions are unsupported by this implementation
    void request1DArray(size_t size) { Log(EError, "request1DArray(): Unsupported!"); }
    void request2DArray(size_t size) { Log(EError, "request2DArray(): Unsupported!"); }
    void setSampleIndex(size_t sampleIndex) { Log(EError, "setSampleIndex(): Unsupported!"); }
    virtual ref<Sampler> clone() { NotImplementedError(); }


    /* ==================================================================== */
    /*                           DRMLT Specifics                            */
    /* ==================================================================== */
    
    /// Set up transition kernels for delayed rejection
    virtual void configureStages() { NotImplementedError(); }

    /// Tell sampler that it is now in reverse mode (EGreen only)
    virtual void setReverse(bool value) { NotImplementedError(); }

    /// Modify sampler behavior in light tracing case (EMMLT only)
    virtual void handleLightTracing() { NotImplementedError(); }

    /// Set all stages to identity map, EMMLT and emitter sampler only
    void setStagesToIdentity();

    /// Fill remaining of current state space to match maximum dimensions after a replay
    inline void fillReplay() {
        while(m_uCurrent.size() < m_maxDim) {
            m_uCurrent.push_back(m_random->nextFloat());
        }
    }
    
    /// Propose new primary samples based on stage (called from primarySample())
    virtual void fillSpace(bool isFirst);

    /// Compute MH transition ratio (only EMira)
    virtual Float getTransitionRatio(Float pLarge) { return 1.f; }

    /// PSS wrap-around helper
    inline Float wrap(Float y) const {
        y = y > 1 ? 2.f - y : (y <= 0 ? std::abs(y) : y);
        SAssert(y >= 0.f && y <= 1.f);
        return y;
    };
    
    /** Increment delayed rejection stage. if pure emitter path and 
     * fixEmitterPath ==true, do light tracing.
    */
    inline void nextStage(bool lightTracing = false) {
        m_sampleIndex = 0;
        m_kernelHierarchy.isFirst = false;
        m_kernelHierarchy.isLightTracing = lightTracing;
    }

    MTS_DECLARE_CLASS()

protected:

	/* ==================================================================== */
    /*                           General                                    */
    /* ==================================================================== */
	
    /// Virtual destructor
    ~DRMLTSampler() override;
    

    ref<Random> m_random;                    // Random number generator
    bool m_replay = false;                   // replay mode
    size_t m_maxDim = 80;                    // maximum dimension TODO: should not be hardcoded
    bool m_largeStep;                        // large step probability

	/* ==================================================================== */
    /*                           DRMLT Specifics                            */
    /* ==================================================================== */

    /**
     * Kernel hierarchy for delayed rejection MLT
    */
    struct KernelHierarchy {
        std::unique_ptr<TransitionKernel> stage1;         // first stage kernel
        std::unique_ptr<TransitionKernel> stage2;         // second stage kernel
        std::unique_ptr<TransitionKernel> stageLT;        // light tracing kernel
        bool isLightTracing = false;                      // is light tracing (in mmlt)
        bool isFirst = true;                              // is first stage
        bool isReverse = false;                           // Is reverse stage, Green & Mira (2001) only
        
        /// get current stage transition kernel
        const TransitionKernel &getCurrent() {
            if (isFirst) return (*stage1);
            else {
                if (isLightTracing) return (*stageLT);
                else return (*stage2);
            }
        }
    };

    DRMLTConfiguration::EType m_type;        // DR algo; green, mira, orbital

    KernelHierarchy m_kernelHierarchy;       // Kernel hierarchy

    const Float m_s1 = 1.f / 1024.f;         // Kelemen - inner bound
    const Float m_s2 = 1.f / 64.f;           // Kelemen - outer bound
    Float m_sigma;                           // Gaussian std dev
    const Float m_rho = std::exp(-0.25f);    // warped cauchy coefficient
    const Float m_kelemenScale = 1.9f;       // EOrbital - compensate for pairwise sampling v.s. i.i.d.
    Float m_scaleSecond;                     // second stage scale

    /* Proposed and current pss states */
    std::pair<std::vector<Float>, std::vector<Float>> m_uProposed;
    std::vector<Float> m_uCurrent;

    size_t m_dimStage1 = 0, m_dimStage2 = 0; // Dimension of proposed states

};

/**
 * Sampler implementation of DRMLT using Green & Mira (2001).
*/
class GreenDRMLTSampler : public DRMLTSampler {
public:
    /// Constructors
    explicit GreenDRMLTSampler(const DRMLTConfiguration &conf) 
        : DRMLTSampler(conf) { configureStages(); }
    explicit GreenDRMLTSampler(DRMLTSampler *sampler) 
        : DRMLTSampler(sampler) {  configureStages(); }

    /// Clone
    ref<Sampler> clone() override;

    /// Set kelemen  as first stage  and small gaussian as second
    void configureStages() override;
    void handleLightTracing() override;
    
    /// Set whether we're computing virtual reverse move from z->y*
    void setReverse(bool value) override {
        m_sampleIndex = 0;
        m_kernelHierarchy.isReverse = value;
    }

    /// Need to account for reverse mode
    Float primarySample(size_t k) override;


    MTS_DECLARE_CLASS()
};


/**
 * Sampler implementation of DRMLT using Tierney & Mira (1999).
*/
class MiraDRMLTSampler : public DRMLTSampler {
public:

    /// Constructors
    explicit MiraDRMLTSampler(const DRMLTConfiguration &conf)
        : DRMLTSampler(conf) { configureStages(); }
    explicit MiraDRMLTSampler(DRMLTSampler *sampler)
        : DRMLTSampler(sampler) { configureStages(); }

    /// Clone
    ref<Sampler> clone() override;

    /// Set kelemen  as first stage  and small gaussian as second
    void configureStages() override;
    void handleLightTracing() override;

    /// Compute problematic transition kernel ratio
    Float getTransitionRatio(Float pLarge) override;


    MTS_DECLARE_CLASS()
};

/**
 * Sampler implementation of DRMLT using Tierney & Mira (1999)formulation
 * in conjunction with Pairwise Orbital mutation.
*/
class OrbitalDRMLTSampler : public DRMLTSampler {
public:

    /// Constructors
    explicit OrbitalDRMLTSampler(const DRMLTConfiguration &conf)
        : DRMLTSampler(conf) { configureStages(); }
    explicit OrbitalDRMLTSampler(DRMLTSampler *sampler)
        : DRMLTSampler(sampler) { configureStages(); }

    /// Clone
    ref<Sampler> clone() override;

    /// Set kelemen  as first stage  and orbital as second
    void configureStages() override;
    void handleLightTracing() override;

    /// Fill states in pairwise fashion
    void fillSpace(bool isFirst) override;
    
    MTS_DECLARE_CLASS()
};

MTS_NAMESPACE_END

#endif /* __DRMLT_SAMPLER_H */