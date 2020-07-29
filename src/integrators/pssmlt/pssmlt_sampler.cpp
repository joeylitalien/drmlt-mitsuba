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

#include "pssmlt_sampler.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                    Constructors, Serialize, Clone, ...               */
/* ==================================================================== */

PSSMLTSampler::PSSMLTSampler(const PSSMLTConfiguration &config) : Sampler(Properties()) {
	m_random = new Random();
	m_s1 = config.mutationSizeLow;
	m_s2 = config.mutationSizeHigh;
    m_sigma = config.sigma;
	configure();
}

PSSMLTSampler::PSSMLTSampler(PSSMLTSampler *sampler) : Sampler(Properties()),
	m_random(sampler->m_random) {
	m_s1 = sampler->m_s1;
	m_s2 = sampler->m_s2;
    m_sigma = sampler->m_sigma;
	configure();
}

PSSMLTSampler::PSSMLTSampler(Stream *stream, InstanceManager *manager)
	: Sampler(stream, manager) {
	m_random = static_cast<Random *>(manager->getInstance(stream));
	m_s1 = stream->readFloat();
	m_s2 = stream->readFloat();
	m_sigma = stream->readFloat();
	configure();
}

PSSMLTSampler::~PSSMLTSampler() { }

void PSSMLTSampler::serialize(Stream *stream, InstanceManager *manager) const {
	Sampler::serialize(stream, manager);
	manager->serialize(stream, m_random.get());
	stream->writeFloat(m_s1);
	stream->writeFloat(m_s2);
	stream->writeFloat(m_sigma);
}

ref<Sampler> PSSMLTSampler::clone() {
	ref<PSSMLTSampler> sampler = new PSSMLTSampler(this);
	sampler->m_sampleCount = m_sampleCount;
	sampler->m_sampleIndex = m_sampleIndex;
	sampler->m_random = new Random(m_random);
	return sampler.get();
}

std::string PSSMLTSampler::toString() const {
	std::ostringstream oss;
	oss << "PSSMLTSampler[" << endl
		<< "  sampleCount = " << m_sampleCount << endl
		<< "]";
	return oss.str();
}
void PSSMLTSampler::configure() {
	m_logRatio = -math::fastlog(m_s2/m_s1);
	m_time = 0;
	m_largeStepTime = 0;
	m_largeStep = false;
	m_sampleIndex = 0;
	m_sampleCount = 0;
}

/* ==================================================================== */
/*                           Sampling Specifics                         */
/* ==================================================================== */

/** 
 * Set the proposed state as new current state, then clear the backup
*/
void PSSMLTSampler::accept() {
	if (m_largeStep)
		m_largeStepTime = m_time;
	m_time++;
	m_backup.clear();
	m_sampleIndex = 0;
}

/** 
 * Reset the current state
*/
void PSSMLTSampler::reset() {
	m_time = m_sampleIndex = m_largeStepTime = 0;
	m_u.clear();
}

/** 
 * Set the state back to the backuped one
*/
void PSSMLTSampler::reject() {
	for (size_t i=0; i<m_backup.size(); ++i)
		m_u[m_backup[i].first] = m_backup[i].second;
	m_backup.clear();
	m_sampleIndex = 0;
}

/** 
 * EMira, EOrbital - Return the k-th component of the proposed state 
 * If k == 0, all components are generated. All further call just return 
 * those precomputed value
*/
Float PSSMLTSampler::primarySample(size_t i) {
    
	/* If replay mode, return from m_random directly */
	if(m_replay) {
        m_u.push_back(SampleStruct(m_random->nextFloat()));
        return m_u[i].value;
    }
	
	/* Fill until maximum dimension is reached */
    if(i == 0) {
        for (auto k = 0; k < m_maxDim; k++) {
			/* Should not be called more often than max dim */
            if(k > m_u.size()) {
                SLog(EError, "Impossible!!!");
            }

			/* if need new dimension generate new rng */
            if(k == m_u.size()) {
                m_u.push_back(SampleStruct(m_random->nextFloat()));
            } 
			/* Otherwise */
			else {
        		/* Large step -> uniform */
                if (m_largeStep) {
                    m_backup.push_back(std::pair<size_t, SampleStruct>(k, m_u[k]));
                    m_u[k].modify = m_time;
                    m_u[k].value = m_random->nextFloat();
                } 
				/* Otherwise, generate from mutate */
				else {
                    m_backup.push_back(std::pair<size_t, SampleStruct>(k, m_u[k]));
                    m_u[k].value = mutate(m_u[k].value);
                    m_u[k].modify++;
                }
            }
        }
    }
    if(i >= m_maxDim) {
        SLog(EError, "out of bound... k=%i", i);
    }

    return m_u[i].value;
}

Float PSSMLTSampler::next1D() {
	return primarySample(m_sampleIndex++);
}

Point2 PSSMLTSampler::next2D() {
	/// Enforce a specific order of evaluation
	Float value1 = primarySample(m_sampleIndex++);
	Float value2 = primarySample(m_sampleIndex++);
	return Point2(value1, value2);
}


MTS_IMPLEMENT_CLASS_S(PSSMLTSampler, false, Sampler)
MTS_NAMESPACE_END
