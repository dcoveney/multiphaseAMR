
#include <AMReX_LevelBld.H>
#include <Euler.H>

using namespace amrex;

class LevelBldAdv
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
				  const DistributionMapping& dm,
                                  Real            time) override;
};

LevelBldAdv Adv_bld;

LevelBld*
getLevelBld ()
{
    return &Adv_bld;
}

void
LevelBldAdv::variableSetUp ()
{
    Euler::variableSetUp();
}

void
LevelBldAdv::variableCleanUp ()
{
    Euler::variableCleanUp();
}

AmrLevel*
LevelBldAdv::operator() ()
{
    return new Euler;
}

AmrLevel*
LevelBldAdv::operator() (Amr&            papa,
	   	         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new Euler(papa, lev, level_geom, ba, dm, time);
}
