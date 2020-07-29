#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/bidir/path.h>

#ifndef MITSUBA_PSSMLT_UTILS_H
#define MITSUBA_PSSMLT_UTILS_H

MTS_NAMESPACE_BEGIN

/**
 * Maximum dimensionality for each sampler (in PSS)
*/
struct MaxDim {
    int sensor;
    int emitter;
    int direct;
};

/**
 * Return the maximum number of PSS components for each 
 * sampler depending on the technique used and scene properties
 * such as: Participating media, rough dielectric and 
 * Russian roulette.
*/
MaxDim findMaxDimensions(const Scene* scene, 
                            int maxDepth, 
                            int rrDepth, 
                            int depth,
                            PathSampler::ETechnique tech, 
                            bool useDirectSampling) {
    
    /* Checking some scene properties */
    bool hasMedium = scene->hasMedia();
    bool hasRoughDielectric = [&]() -> bool {
        auto shapes = scene->getShapes();
        auto foundRoughDielectric = false;
        for (auto &s: shapes) {
            if (s->getBSDF() != nullptr) {
                auto bsdf = s->getBSDF();
                foundRoughDielectric |= bsdf->getClass()->getName() == "RoughDielectric";
            }
        }
        return foundRoughDielectric;
    }();
    bool hasRR = rrDepth < maxDepth;

    /* Compute the offset per dimension that we need to apply */
    int offsetMedium = hasMedium ? 1 : 0;
    int offsetRoughDielectric = hasRoughDielectric ? 1 : 0;
    int offsetRR = hasRR ? 1 : 0;
    

    /* return the right amount depending of the technique */
    if (tech == PathSampler::EMMLT) {
        int maxDim = (depth + 2) * (3); 
        if (maxDim % 2 == 1) {
            maxDim++;
        }
        return MaxDim{ maxDim, maxDim, 1};
    } else if (tech == PathSampler::EUnidirectional) {
        SAssert(maxDepth != -1);
        int maxDim = (maxDepth + 2) * (4 + offsetMedium + offsetRoughDielectric + offsetRR);
        if (maxDim % 2 == 1) {
            maxDim++;
        }
        return MaxDim{maxDim, 0, 0};
    } else if (tech == PathSampler::EBidirectional) {
        SAssert(maxDepth != -1);
        int maxDim = (maxDepth + 2) * (2 + offsetMedium + offsetRoughDielectric + offsetRR);
        if (maxDim % 2 == 1) {
            maxDim++;
        }
        return MaxDim{maxDim, maxDim, useDirectSampling ? maxDepth : 0};
    }
}

MTS_NAMESPACE_END

#endif //MITSUBA_PSSMLT_UTILS_H