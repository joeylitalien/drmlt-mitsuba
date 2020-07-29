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


#if !defined(__PSSMLT_SAMPLER_H)
#define __PSSMLT_SAMPLER_H

#include <mitsuba/render/sampler.h>
#include <mitsuba/core/random.h>
#include "pssmlt.h"

MTS_NAMESPACE_BEGIN

/**
 * Sampler implementation as described in
 * 'A Simple and Robust Mutation Strategy for the
 * Metropolis Light Transport Algorithm' by Kelemen et al.
*/
class PSSMLTSampler : public Sampler {
public:

	/* ==================================================================== */
    /*                           General                                    */
    /* ==================================================================== */

	/// Construct a new MLT sampler
	PSSMLTSampler(const PSSMLTConfiguration &conf);

	/* Construct a new sampler, which operates on the
	same random number generator as a sampler.*/
	PSSMLTSampler(PSSMLTSampler *sampler);

	/// Unserialize from a binary data stream
	PSSMLTSampler(Stream *stream, InstanceManager *manager);

	/// Set up the internal state
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

    /// Set replay mode
    inline void setReplay(bool value) { m_replay = value; }

    /// Set maximum dimension of primary sample vector
    inline void setMaxDim(size_t dim) { m_maxDim = dim; }

	/// Return a primary sample
	Float primarySample(size_t i);

	/// Reset (& start with a large mutation)
	void reset();

	/// Accept a mutation
	void accept();

	/// Reject a mutation
	void reject();

	/// Replace the underlying random number generator
	inline void setRandom(Random *random) { m_random = random; }

	/// Return the underlying random number generator
	inline Random *getRandom() { return m_random; }

	/* The following functions do nothing in this implementation */
	virtual void advance() { }
	virtual void generate(const Point2i &pos) { }

	/* The following functions are unsupported by this implementation */
	void request1DArray(size_t size) { Log(EError, "request1DArray(): Unsupported!"); }
	void request2DArray(size_t size) { Log(EError, "request2DArray(): Unsupported!"); }
	void setSampleIndex(size_t sampleIndex) { Log(EError, "setSampleIndex(): Unsupported!"); }
	ref<Sampler> clone();

	/* ==================================================================== */
    /*                           PSSMLT Specifics                           */
    /* ==================================================================== */


	/// Set whether to use Kelemen mutation or not; if false, use Gaussian
	inline void setMutationType(bool value) {
	    m_useKelemen = value;
	}

	/// 1D mutation routine
	inline Float mutate(Float value) {
	    if (m_useKelemen) {
            Float sample = m_random->nextFloat();
            bool add;

            if (sample < 0.5f) {
                add = true;
                sample *= 2.0f;
            } else {
                add = false;
                sample = 2.0f * (sample - 0.5f);
            }

            Float dv = m_s2 * math::fastexp(sample * m_logRatio);
            if (add) {
                value += dv;
                if (value > 1)
                    value -= 1;
            } else {
                value -= dv;
                if (value < 0)
                    value += 1;
            }
        } else {
            Float tmp1 = std::sqrt(-2 * std::log(1 - m_random->nextFloat()));
            Float dv = tmp1 * std::cos(2 * M_PI * m_random->nextFloat());
            value = math::modulo(value + m_sigma * dv, 1.0f);
        }

		return value;
	}


	MTS_DECLARE_CLASS()

protected:

	/* ==================================================================== */
    /*                           General                                    */
    /* ==================================================================== */

	/// Virtual destructor
	virtual ~PSSMLTSampler();

    ref<Random> m_random;               // Random number generator
    bool m_replay = false;              // replay mode
    size_t m_maxDim = 80;               // maximum dimension TODO: should not be hardcoded
    bool m_largeStep;                   // large step probability

	/* ==================================================================== */
    /*                           PSSMLT Specifics                           */
    /* ==================================================================== */

	/**
	 * keep count of last use of a pss dimension
	*/
	struct SampleStruct {
		Float value;
		size_t modify;
		inline SampleStruct(Float value) : value(value), modify(0) { }
	};

	bool m_useKelemen;                                      // Use Kelemen of Gaussian
	Float m_s1, m_s2, m_logRatio;                           // Kelemen mutation related quantities
	Float m_sigma;                                          // Gaussian mutation std dev

	std::vector<SampleStruct> m_u;	                        // state
	std::vector<std::pair<size_t, SampleStruct> > m_backup;	// state backup
	size_t m_time, m_largeStepTime;                         // Kelemen last large step

};

MTS_NAMESPACE_END

#endif /* __PSSMLT_SAMPLER_H */
