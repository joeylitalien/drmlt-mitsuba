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

#include "drmlt_sampler.h"

#include <memory>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                    Constructors, Serialize, Clone, ...               */
/* ==================================================================== */

DRMLTSampler::DRMLTSampler(const DRMLTConfiguration &config) : Sampler(Properties()) {
    m_random = new Random();
    m_type = config.type;
    m_sigma = config.sigma;
    m_scaleSecond = config.scaleSecond;
    configure();
}

DRMLTSampler::DRMLTSampler(DRMLTSampler *sampler) : Sampler(Properties()),
                                                    m_random(sampler->m_random) {
    m_sigma = sampler->m_sigma;
    m_type = sampler->m_type;
    m_scaleSecond = sampler->m_scaleSecond;
    configure();
}

DRMLTSampler::DRMLTSampler(Stream *stream, InstanceManager *manager)
    : Sampler(stream, manager) {
    m_random = static_cast<Random *>(manager->getInstance(stream));
    m_type = static_cast<DRMLTConfiguration::EType>(stream->readShort());
    m_sigma = stream->readFloat();
    m_scaleSecond = stream->readFloat();
    configure();
}

DRMLTSampler::~DRMLTSampler() {}

void DRMLTSampler::serialize(Stream *stream, InstanceManager *manager) const {
    Sampler::serialize(stream, manager);
    manager->serialize(stream, m_random.get());
    stream->writeShort(m_type);
    stream->writeFloat(m_sigma);
    stream->writeFloat(m_scaleSecond);
}

ref<Sampler> GreenDRMLTSampler::clone() {
    ref<GreenDRMLTSampler> sampler = new GreenDRMLTSampler(this);
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_sampleIndex = m_sampleIndex;
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

ref<Sampler> MiraDRMLTSampler::clone() {
    ref<MiraDRMLTSampler> sampler = new MiraDRMLTSampler(this);
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_sampleIndex = m_sampleIndex;
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

ref<Sampler> OrbitalDRMLTSampler::clone() {
    ref<OrbitalDRMLTSampler> sampler = new OrbitalDRMLTSampler(this);
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_sampleIndex = m_sampleIndex;
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

std::string DRMLTSampler::toString() const {
    std::ostringstream oss;
    oss << "DRMLTSampler[" << endl
        << "  sampleCount = " << m_sampleCount << endl
        << "]";
    return oss.str();
}

void DRMLTSampler::configure() {
    m_largeStep = false;
    m_sampleIndex = 0;
    m_sampleCount = 0;
    m_kernelHierarchy.isFirst = true;
    m_kernelHierarchy.isReverse = false;
}


/* ==================================================================== */
/*                           Stages Specifics                           */
/* ==================================================================== */

/**
 * Set all stage to identity. (EMMLT - direct sampler) 
*/
void DRMLTSampler::setStagesToIdentity() { 
    m_kernelHierarchy.stage1  = std::make_unique<IdentityKernel>();
    m_kernelHierarchy.stage2  = std::make_unique<IdentityKernel>();
    m_kernelHierarchy.stageLT = std::make_unique<IdentityKernel>();
}

/**
 * EGreen - Set first stage to Kelemen and second to scaled Gaussian
*/
void GreenDRMLTSampler::configureStages() {
    m_kernelHierarchy.stage1 = std::make_unique<KelemenKernel>(m_s1, m_s2);
    m_kernelHierarchy.stage2 = std::make_unique<GaussianKernel>(m_scaleSecond * m_sigma);
}

/**
 * EGreen - if fixed emitter (EMMLT), change second stage to identity in emitter sampler,  
 * otherwise same as configureStages
*/
void GreenDRMLTSampler::handleLightTracing() {
    m_kernelHierarchy.stage1  = std::make_unique<KelemenKernel>(m_s1, m_s2);
    m_kernelHierarchy.stage2  = std::make_unique<IdentityKernel>();
    m_kernelHierarchy.stageLT = std::make_unique<GaussianKernel>(m_scaleSecond * m_sigma);
}

/**
 * EMira - Set first stage to Kelemen and second to scaled Gaussian 
*/
void MiraDRMLTSampler::configureStages() {
    m_kernelHierarchy.stage1 = std::make_unique<KelemenKernel>(m_s1, m_s2);
    m_kernelHierarchy.stage2 = std::make_unique<GaussianKernel>(m_scaleSecond * m_sigma);
}

/**
 * EMira - if fixed emitter (EMMLT), change second stage to identity in emitter sampler. 
 * otherwise same as configureStages
*/
void MiraDRMLTSampler::handleLightTracing() {
    m_kernelHierarchy.stage1 = std::make_unique<KelemenKernel>(m_s1, m_s2);
    m_kernelHierarchy.stage2  = std::make_unique<IdentityKernel>();
    m_kernelHierarchy.stageLT = std::make_unique<GaussianKernel>(m_scaleSecond * m_sigma);
}

/**
 * EOrbital - Set first stage to Kelemen, scaled for pairwise sampling, 
 * and second a warped cauchy
*/
void OrbitalDRMLTSampler::configureStages() {
     // compensate pairwise sampling against product of i.i.d.
    Float s1 = m_s1 * m_kelemenScale;
    Float s2 = m_s2 * m_kelemenScale;
    m_kernelHierarchy.stage1 = std::make_unique<KelemenKernel>(s1, s2);
    m_kernelHierarchy.stage2 = std::make_unique<WrappedCauchyKernel>(m_rho);
}

/**
 * EOrbital - if fixed emitter (EMMLT), change second stage to identity in emitter sampler. 
 * otherwise same as configureStages 
*/
void OrbitalDRMLTSampler::handleLightTracing() {
    // compensate pairwise sampling against product of i.i.d.
    Float s1 = m_s1 * m_kelemenScale;
    Float s2 = m_s2 * m_kelemenScale;
    m_kernelHierarchy.stage1  = std::make_unique<KelemenKernel>(s1, s2);
    m_kernelHierarchy.stage2  = std::make_unique<IdentityKernel>();
    m_kernelHierarchy.stageLT = std::make_unique<WrappedCauchyKernel>(m_rho);
}


/* ==================================================================== */
/*                           Sampling Specifics                         */
/* ==================================================================== */


/** 
 * Set the accepted proposed state as new current state, then reset to default values 
 * and empty proposed states
*/
void DRMLTSampler::accept(bool acceptFirst) {
    m_uCurrent = acceptFirst ? m_uProposed.first : m_uProposed.second;
    for (Float & i : m_uCurrent) i = wrap(i);

    /* And "reset" (same as reject) */
    m_sampleIndex = 0;
    m_kernelHierarchy.isFirst = true;
    m_kernelHierarchy.isReverse = false;
    m_uProposed.first.clear(); m_uProposed.second.clear();
    m_dimStage1 = 0; m_dimStage2 = 0;
}

/** 
 * Reset to default value, empty all states
*/
void DRMLTSampler::reset() {
    m_sampleIndex = 0;
    m_kernelHierarchy.isFirst = true;
    m_kernelHierarchy.isReverse = false;
    m_uCurrent.clear();
    m_uProposed.first.clear(); m_uProposed.second.clear();
    m_dimStage1 = 0; m_dimStage2 = 0;
}

/** 
 * Reset to default values and empty proposed states
*/
void DRMLTSampler::reject() {
    m_sampleIndex = 0;
    m_kernelHierarchy.isFirst = true;
    m_kernelHierarchy.isReverse = false;
    m_uProposed.first.clear(); m_uProposed.second.clear();
    m_dimStage1 = 0; m_dimStage2 = 0;
}


/** 
 * EMira, EOrbital - Return the k-th component of the proposed state 
 * If k == 0, all components are generated and saved in the appropriated 
 * m_uProposed, depending on the stage. All further call just return from the 
 * approriate m_uProposed.
*/
Float DRMLTSampler::primarySample(size_t k) {
    auto &uProposed = m_kernelHierarchy.isFirst ?
                        m_uProposed.first :
                        m_uProposed.second;

    
    m_kernelHierarchy.isFirst ? m_dimStage1 = std::max(k, m_dimStage1)
                                : m_dimStage2 = std::max(k, m_dimStage2);

    /* If first call, clear proposed state to fill it later. */
    if (k == 0) {
        uProposed.clear();
    }

    /* If replay mode, return from m_random directly */
    if (m_replay) {
        uProposed.push_back(m_random->nextFloat());
        return wrap(uProposed[k]);
    }

    /* Otherwise, call fillSpace to fill the propose state */
    if (k == 0) {
        fillSpace(m_kernelHierarchy.isFirst);
    }

    /* Should not exceed m_maxDim */
    if (k > m_maxDim) {
        SLog(EError, "Exceeded maximum dimension: k = %i (Max = %i)", k, m_maxDim);
    }

    /* Return the wrapped appropriate component */
    return wrap(uProposed[k]);
}

/** 
 * EGreen - Return the k-th component of the proposed state.
 * If we are in reversed mode, return y^* = z-(y-x).
 * Otherwise, same as the base function. See DRMLT (2020), 
 * Sec 4.4
*/
Float GreenDRMLTSampler::primarySample(size_t k) {
    auto &uProposed = m_kernelHierarchy.isFirst ? m_uProposed.first
                                                : m_uProposed.second;

    m_kernelHierarchy.isFirst ? m_dimStage1 = std::max(k, m_dimStage1)
                                : m_dimStage2 = std::max(k, m_dimStage2);

    /* If not in reverse mode and this is the first call, 
    clear proposed state to fill it later. */
    if (k == 0 && !m_kernelHierarchy.isReverse) {
        uProposed.clear();
    }

    /* If replay mode, return from m_random directly */
    if (m_replay) {
        uProposed.push_back(m_random->nextFloat());
        return wrap(uProposed[k]);
    }

    /* If in reverse mode return y^* = z-(y-x). See DRMLT Sec. 4.4  */
    if (m_kernelHierarchy.isReverse) {
        auto du =  m_uProposed.first[k] - m_uCurrent[k];
        return wrap(m_uProposed.second[k] - du);
    } 
    /* Otherwise, call fillSpace to fill the propose state */
    else if (k == 0) {
        fillSpace(m_kernelHierarchy.isFirst);
    }

    /* Should not exceed m_maxDim */
    if (k > m_maxDim) {
        SLog(EError, "Exceeded maximum dimension: k = %i (Max = %i)", k, m_maxDim);
    }

    /* And return the appropriate component */
    return wrap(uProposed[k]);
}

/** 
 * EGreen, EMira - fill the appropriate proposed state in a 
 * i.i.d. fashion
*/
void DRMLTSampler::fillSpace(bool isFirst) {
    auto &uProposed = isFirst ? m_uProposed.first : m_uProposed.second;

    /* Fill until maximum dimension is reached */
    for (auto i = 0; i < m_maxDim; i++) {
        /* Large step -> uniform */
        if (m_largeStep) {
            SAssert(isFirst);
            uProposed.push_back(m_random->nextFloat());
        /* Identity -> don't move (fixed emmiter) */
        } else if (m_kernelHierarchy.getCurrent().isIdentity()) {
            uProposed.push_back(m_uCurrent[i]);
        }
        /* Otherwise, generate i-th component of proposal state from the right transition kernel */
        else {
            SAssert(!m_uCurrent.empty());
            uProposed.push_back(m_uCurrent[i] + m_kernelHierarchy.getCurrent().sample(m_random));
        }
    }
}

/** 
 * EOrbital - fill the appropriate proposed state in a pairwise fashion
 * If we are in the second stage construct samples according to the 
 * orbital mutation framework. See DRMLT (2020), Sec 4.3 
*/
void OrbitalDRMLTSampler::fillSpace(bool isFirst) {
    auto &uProposed = isFirst ? m_uProposed.first : m_uProposed.second;

     /* Fill until maximum dimension is reached */
    for (auto i = 0; i < m_maxDim; i++) {
        /* Largestep -> uniform */
        if (m_largeStep) {
            SAssert(isFirst);
            uProposed.push_back(m_random->nextFloat());
        } 
        /* Identity -> don't move (fixed emmiter) */
        else if (m_kernelHierarchy.getCurrent().isIdentity()) {
            uProposed.push_back(m_uCurrent[i]);
        } 
        /* First stage - generate pair of samples from a 2D kelemen kernel */
        else if (isFirst) {
            SAssert(m_uCurrent.size() > i + 1);
            auto d = m_kernelHierarchy.getCurrent().sample(m_random);
            auto a = m_random->nextFloat() * 2.f * M_PI;
            uProposed.push_back(m_uCurrent[i] + d * std::cos(a));
            uProposed.push_back(m_uCurrent[i + 1] + d * std::sin(a));
            i++;
        } 
        /* Second stage - generate pair of samples from the orbital mutation */
        else {
            SAssert(m_uCurrent.size() > i + 1);
            SAssert(m_uProposed.first.size() > i + 1);

            /* Main idea:  sample a point near x that is also on an hypersphere of radius ||x-y|| around y
            that way, the q_1(y|z)/q_1(y|x) term diseapear in the acceptance ratio.
            To do so: generate an angular perturbation around the x-y axis with the use of a
            wrapped normal around the hypersphere. we do so coordinate wise and use hyperspherical
            coordinate to do so to get a normalize vector at the end. See DRMLT (2020), Sec 4.4 */

            /* Do two by two... */
            Float theta_i = m_kernelHierarchy.getCurrent().sample(m_random);

            auto du1 = m_uProposed.first[i] - m_uCurrent[i];
            auto du2 = m_uProposed.first[i + 1] - m_uCurrent[i + 1];

            /* Spherical coordinate between the i and i+1 coordinate of x-y = -du1 */
            Float norm = sqrt(du1 * du1 + du2 * du2);
            Float mu_i = math::safe_acos(-du1 / norm);
            if (-du2 < 0)
                mu_i = 2.f * M_PI - mu_i;

            /*  Add projection to main axis */
            auto c1 = m_uProposed.first[i] + cos(theta_i + mu_i) * norm;
            auto c2 = m_uProposed.first[i + 1] + sin(theta_i + mu_i) * norm;

            uProposed.push_back(c1);
            uProposed.push_back(c2);
            i++;
        }
    }
}


/**
 *  EMira - Compute Q1(z | y) / Q1(x | y) 
*/
Float MiraDRMLTSampler::getTransitionRatio(Float pLarge) {
    if (m_kernelHierarchy.stage1->isIdentity()) {
        return 1.f;
    }

    auto dimStage = std::max(m_dimStage1, m_dimStage2);
    Float num = 0.f, denum = 0.f;
    /* Log for numerical stability */
    for (size_t i = 0; i < dimStage; i++) {
        num   += m_kernelHierarchy.stage1->logPdf(m_uProposed.second[i] - m_uProposed.first[i]);
        denum += m_kernelHierarchy.stage1->logPdf(m_uCurrent[i] - m_uProposed.first[i]);
    }

    return std::exp(num - denum);
}

Float DRMLTSampler::next1D() {
    return primarySample(m_sampleIndex++);
}

Point2 DRMLTSampler::next2D() {
    /* Enforce a specific order of evaluation*/
    Float value1 = primarySample(m_sampleIndex++);
    Float value2 = primarySample(m_sampleIndex++);
    return Point2(value1, value2);
}



MTS_IMPLEMENT_CLASS_S(DRMLTSampler, true, Sampler)
MTS_IMPLEMENT_CLASS(GreenDRMLTSampler, false, DRMLTSampler)
MTS_IMPLEMENT_CLASS(MiraDRMLTSampler, false, DRMLTSampler)
MTS_IMPLEMENT_CLASS(OrbitalDRMLTSampler, false, DRMLTSampler)
MTS_NAMESPACE_END