#include <memory>

#if !defined(__DRMLT_TRANSITION_H)
#define __DRMLT_TRANSITION_H

#define M_SQRT1_2PI 0.39894228040143267794

MTS_NAMESPACE_BEGIN

/**
 *  Specifies the type of Markov transition kernel
*/
enum KernelType {
    EGaussian,      // 1D Gaussian 
    EKelemen,       // 1D Kelemen (hole)
    EIdentity,      // 1D Identity (dirac)
    EWrappedCauchy  // Circular - Wrapped Cauchy
};

/**
 * \brief Abstract class for Markov transition kernel
*/
class TransitionKernel {
public:
    /// Constructor
    explicit TransitionKernel(const KernelType _type) : type(_type) {}
    
    /// Destructor
    virtual ~TransitionKernel() = default;
    
    /// Generate a random sample
    virtual Float sample(Random *random) const = 0;
    
    /// Evaluate the pdf at du
    virtual Float pdf(Float du) const = 0;
    
    /// Evaluate the log pdf at du
    virtual Float logPdf(Float du) const = 0;
    
    /// Return the type of kernel
    KernelType getType() const { return type; };
    
    /// Return true if the type is identity
    virtual bool isIdentity() const { return false; }

private:
    const KernelType type;  // type of the kernel
};


/**
 * \brief Gaussian transition kernel with zero-mean
*/
class GaussianKernel : public TransitionKernel {
public:
    explicit GaussianKernel(const Float _sigma)
        : TransitionKernel(KernelType::EGaussian),
        sigma(_sigma), sigma_log(std::log(_sigma*_sigma)) {
    }

    inline Float sample(Random *random) const override {
        /// Box-Muller transform
        Float tmp = std::sqrt(-2.0f * std::log(1 - random->nextFloat()));
        Float dv = tmp * std::cos(2.0f * M_PI * random->nextFloat());
        return dv * sigma;
    }

    inline Float pdf(Float du) const override {
        Float invSigma = 1.0f / sigma;
        Float d = du;
        Float p = M_SQRT1_2PI * invSigma * math::fastexp(-0.5f * d * d * invSigma * invSigma);
        return p;
    }

    inline Float logPdf(Float du) const override {
        Float invSigma = 1.0f / sigma;
        Float r = invSigma * du;
        return -0.5 * ( r*r + std::log(2*M_PI) + sigma_log);
    }

private:
    const Float sigma;      // std dev
    const Float sigma_log;  // log variance
};


/**
 * \brief Kelemen-style transition kernel
*/
class KelemenKernel : public TransitionKernel {
public:
    KelemenKernel(const Float _s1, const Float _s2)
        : TransitionKernel(KernelType::EKelemen), s1(_s1), s2(_s2) {
        logRatio = -math::fastlog(_s2 / _s1);
    }

    inline Float sample(Random *random) const override {
        Float sample = random->nextFloat();
        int sign;

        if (sample < 0.5f) {
            sign = 1;
            sample *= 2.0f;
        } else {
            sign = -1;
            sample = 2.0f * (sample - 0.5f);
        }

        Float dv = s2 * math::fastexp((1-sample) * logRatio);
        return dv * sign;
    }

    inline Float pdf(Float du) const override {
        Float d = std::abs(du);
        if (d < s1 || d > s2) return 0.0f;
        Float p = 1.0f / (2.0f * d * (-logRatio));
        return p;
    }

    inline Float logPdf(Float du) const override {
        return math::fastlog(pdf(du));
    }

private:
    const Float s1, s2;  // Kelemen parameters
    Float logRatio;      // -log(s2/s1)
};


/**
 * \brief Identity transition kernel
*/
class IdentityKernel : public TransitionKernel {
public:
    IdentityKernel() : TransitionKernel(KernelType::EIdentity) {}
    inline Float sample(Random *random) const override { return 0.0f; }
    inline Float pdf(Float du) const override { return 1.0f; }
    inline Float logPdf(Float du) const override { return 0.0f; }
    
    /* Is the identity */
    bool isIdentity() const override { return true; }
};


/**
 * \brief Wrapped cauchy transition kernel with zero-mean
 * See Statistical Analysis of Circular Data: Fisher (1993)
 * and DRMLT (2020) Sec. 4.3
*/
class WrappedCauchyKernel : public TransitionKernel {
public:
    explicit WrappedCauchyKernel(const Float _rho)
        : TransitionKernel(KernelType::EWrappedCauchy),
        rho(_rho), dispersion(2.0*_rho/(1.0+_rho*_rho)) {
    }

    inline Float sample(Random *random) const override {
        Float sample = random->nextFloat();
        int sign = 1;
        if (sample < 0.5f) {
            sign = 1;
            sample *= 2.0f;
        } else {
            sign = -1;
            sample = 2.0f * (sample - 0.5f);
        }
        /* Inversion of analytic CDF */
        Float V = std::cos(2.0*M_PI*sample);
        Float angle = (V + dispersion)/(1.0 + dispersion*V);
        Float dv = sign*math::safe_acos(angle);

        return dv;
    }

    inline Float pdf(Float du) const override {
        /* DRMLT (2020) Eq. 10 */
        Float rho2 = rho*rho;
        Float d = du;
        Float p = 0.5*M_1_PI *(1.0-rho2)/(1.0+rho2-2.0*rho*cos(d));
        return p;
    }

    inline Float logPdf(Float du) const override {
        return math::fastlog(pdf(du));
    }

private:
    const Float rho;          // "variance" of the circular distribution
    const Float dispersion;   // 2*_rho/(1+_rho^2)
};

MTS_NAMESPACE_END

#endif /* __DRMLT_TRANSITION_H */