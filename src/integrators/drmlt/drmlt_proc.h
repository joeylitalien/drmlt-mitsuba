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

#if !defined(__DRMLT_PROC_H)
#define __DRMLT_PROC_H

#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/bidir/rsampler.h>
#include "drmlt.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

class DRMLTProcess : public ParallelProcess {
public:
    DRMLTProcess(const RenderJob *parent, RenderQueue *queue,
                    const DRMLTConfiguration &config, const Bitmap *directImage,
                    const std::vector<PathSeed> &seeds, ref_vector<ReplayableSampler> rplSamplers);

    void develop();

    /* ParallelProcess impl. */
    void processResult(const WorkResult *wr, bool cancelled);
    ref<WorkProcessor> createWorkProcessor() const;
    void bindResource(const std::string &name, int id);
    EStatus generateWork(WorkUnit *unit, int worker);

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~DRMLTProcess() {}
private:
    ref<const RenderJob> m_job;
    RenderQueue *m_queue;
    const DRMLTConfiguration &m_config;
    const Bitmap *m_directImage;
    ref<Bitmap> m_developBuffer;
    ImageBlock *m_accum;
    ProgressReporter *m_progress;
    const std::vector<PathSeed> &m_seeds;
    ref_vector<ReplayableSampler> m_rplSamplers;
    ref<Mutex> m_resultMutex;
    ref<Film> m_film;
    int m_resultCounter, m_workCounter;
    unsigned int m_refreshTimeout;
    ref<Timer> m_timeoutTimer, m_refreshTimer;
};

MTS_NAMESPACE_END

#endif /* __DRMLT_PROC */
