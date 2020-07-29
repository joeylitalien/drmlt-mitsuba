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

#if !defined(__DRMLT_H)
#define __DRMLT_H

#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters used by
 * the MLT rendering implementation
*/
struct DRMLTConfiguration {

    /* ==================================================================== */
    /*                           General                                    */
    /* ==================================================================== */

    PathSampler::ETechnique technique;  // Path sampling technique; bidirectional, multiplexed or unidirectional
    int maxDepth;
    int rrDepth;
    bool directSampling;
    int directSamples;
    bool separateDirect;
    Float luminance;                   // internal
    int luminanceSamples;
    Float pLarge;
    int workUnits;
    size_t nMutations;                 // internal
    bool kelemenStyleWeights;
    bool twoStage;
    bool firstStage;                   // internal
    int firstStageSizeReduction;
    size_t timeout;
    ref<Bitmap> importanceMap;         // internal
    float averageLuminance;
    bool lightImage;   

    /* ==================================================================== */
    /*                           DRMLT Specifics                            */
    /* ==================================================================== */
    
    /**
     * Specifies the sampling algorithm that is internally used
    */
    enum EType {
        EGreen,     // uses Green & Mira (2001) algorithm
        EMira,      // uses Tierney & Mira (1999) algorithm
        EOrbital,   // uses uses Tierney & Mira (1999) algorithm with Orbital second stage mutation
    };
    
    EType type;                        // DR framework to use: green, mira or orbital
    bool acceptanceMap;                // Generate acceptance map
    bool timidAfterLarge;              // perform a second stage after a large step
    bool fixEmitterPath;               // fix emitter path in second stage
    bool useMixture;                   // use mixture of delayed rejection
    Float sigma;                       // Gaussian Kernel standard deviation
    Float scaleSecond;                 // variance scaling at second stage (scaleSecond*sigma)
    
    
    inline DRMLTConfiguration() {}

    void dump() const {
        SLog(EDebug, "General configuration:");
        SLog(EDebug, "   Path sampling technique  : %s",
            (technique == PathSampler::EBidirectional) ? "bdpt" : 
                ((technique == PathSampler::EUnidirectional) ? "path":  
                    ((technique == PathSampler::EMMLT) ? "mmlt": "wrong technique")));
        SLog(EDebug, "   Maximum path depth          : %i", maxDepth);
        SLog(EDebug, "   Russian roulette depth      : %i", rrDepth);
        SLog(EDebug, "   Direct sampling strategies  : %s", directSampling ? "yes" : "no");
        SLog(EDebug, "   Separate direct illum.      : %s", separateDirect ? "yes" : "no");
        SLog(EDebug, "   Overall MLT image luminance : %f (%i samples)", luminance, luminanceSamples);
        SLog(EDebug, "   Large step probability      : %f", pLarge);
        SLog(EDebug, "   Total number of work units  : %i", workUnits);
        SLog(EDebug, "   Direct illum. samples       : %i", directSamples);
        SLog(EDebug, "   Mutations per work unit     : " SIZE_T_FMT, nMutations);
        SLog(EDebug, "   Kelemen et al. weights      : %s", kelemenStyleWeights ? "yes" : "no");
        SLog(EDebug, "   Two-stage MLT               : %s", twoStage ? "yes" : "no");
        if (twoStage)
            SLog(EDebug, "   First-stage size reduction  : %i", firstStageSizeReduction);
        if (timeout)
            SLog(EDebug, "   Timeout                     : ", SIZE_T_FMT, timeout);
        SLog(EDebug, "   Normalization factor        : %f", averageLuminance);
        SLog(EDebug, "     : %s", lightImage ? "yes" : "no");

        SLog(EDebug, "DRMLT specifics configuration:");
        SLog(EDebug, "   Delayed rejection type         : %s",
            (type == EGreen) ? "green" : 
                ((type == EMira) ? "mira":  
                    ((type == EOrbital) ? "orbital": "wrong type")));
        SLog(EDebug, "   Output acceptance map       : %s", acceptanceMap ? "yes" : "no");
        SLog(EDebug, "   Timid after large           : %s", timidAfterLarge ? "yes" : "no");
        SLog(EDebug, "   Fix emitter subpath         : %s", fixEmitterPath ? "yes" : "no");
        SLog(EDebug, "   Use mixture instead of DR   : %s", useMixture ? "yes" : "no");
        SLog(EDebug, "   Gaussian kernel std         : %f", sigma);
        SLog(EDebug, "   Second stage scale          : %f", scaleSecond);
    }

    inline DRMLTConfiguration(Stream *stream) {
        technique = (PathSampler::ETechnique) stream->readUInt();
        maxDepth = stream->readInt();
        rrDepth = stream->readInt();
        directSampling = stream->readBool();
        directSamples = stream->readInt();
        separateDirect = stream->readBool();
        luminance = stream->readFloat();
        luminanceSamples = stream->readInt();
        pLarge = stream->readFloat();
        workUnits = stream->readInt();
        nMutations = stream->readSize();
        kelemenStyleWeights = stream->readBool();
        twoStage = stream->readBool();
        firstStage = stream->readBool();
        firstStageSizeReduction = stream->readInt();
        timeout = stream->readSize();
        Vector2i size(stream);
        if (size != Vector2i(0)) {
            importanceMap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat, size);
            stream->readFloatArray(importanceMap->getFloatData(),
                                   (size_t) size.x * (size_t) size.y);
        }
        averageLuminance = stream->readFloat();
        lightImage = stream->readBool();
        type = (EType) stream->readUInt();
        acceptanceMap = stream->readBool();
        timidAfterLarge = stream->readBool();
        fixEmitterPath = stream->readBool();
        useMixture = stream->readBool();
        sigma = stream->readFloat();
        scaleSecond = stream->readFloat();
    }

// TODO: update those
    inline void serialize(Stream *stream) const {
        stream->writeUInt((uint32_t) technique);
        stream->writeInt(maxDepth);
        stream->writeInt(rrDepth);
        stream->writeBool(directSampling);
        stream->writeInt(directSamples);
        stream->writeBool(separateDirect);
        stream->writeFloat(luminance);
        stream->writeInt(luminanceSamples);
        stream->writeFloat(pLarge);
        stream->writeInt(workUnits);
        stream->writeSize(nMutations);
        stream->writeBool(kelemenStyleWeights);
        stream->writeBool(twoStage);
        stream->writeBool(firstStage);
        stream->writeInt(firstStageSizeReduction);
        stream->writeSize(timeout);
        if (importanceMap.get()) {
            importanceMap->getSize().serialize(stream);
            stream->writeFloatArray(importanceMap->getFloatData(),
                                    (size_t) importanceMap->getWidth() * (size_t) importanceMap->getHeight());
        } else {
            Vector2i(0, 0).serialize(stream);
        }
        stream->writeFloat(averageLuminance);
        stream->writeBool(lightImage);
        stream->writeUInt((uint32_t) type);
        stream->writeBool(acceptanceMap);
        stream->writeBool(timidAfterLarge);
        stream->writeBool(fixEmitterPath);
        stream->writeBool(useMixture);
        stream->writeFloat(sigma);
        stream->writeFloat(scaleSecond);
    }
};


MTS_NAMESPACE_END

#endif /* __DRMLT_H */
