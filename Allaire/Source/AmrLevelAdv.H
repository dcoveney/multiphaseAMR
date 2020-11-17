#ifndef _MultiphaseAMR_H_
#define _MultiphaseAMR_H_

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX.H>
// #include <.H>
#include <memory>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

enum StateType {Cons_Type = 0, Prim_Type};

//
// AmrLevel-derived class for hyperbolic conservation equations
//

class MultiphaseAMR
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    MultiphaseAMR ();
    //
    //The basic constructor.
    //
    MultiphaseAMR (amrex::Amr&     papa,
	         int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~MultiphaseAMR () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    virtual void checkPoint (const std::string& dir,
			     std::ostream&      os,
			     amrex::VisMF::How  how = amrex::VisMF::NFiles,
			     bool               dump_old = true) override;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () override;
    //
    //Initialize data on this level from another MultiphaseAMR (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    void initialState(amrex::FArrayBox& statein, const amrex::Box& bx,
        amrex::RealBox gridloc, const amrex::Real* dx, int testNumber);


    std::vector<amrex_real> fluxHLLCcell(std::vector<amrex_real>& U_R, std::vector<amrex_real>& U_L,
        amrex::Real gam, int Dim, int Nv);

    std::vector<amrex_real> eqnFluxes(std::vector<amrex_real> U, amrex::Real gam, int Dim, int Nv);

    void MUSCLHancock(const amrex::FArrayBox& statein, const amrex::Box& bx,
    amrex::FArrayBox& ULbar, amrex::FArrayBox& URbar, amrex::FArrayBox& flux,
            const amrex::Real* dx, double dt, amrex::Real gam, int Dim);
    void primMUSCLThirdOrder(const amrex::FArrayBox& statein, const amrex::Box& bx,
    amrex::FArrayBox& ULbar, amrex::FArrayBox& URbar, amrex::FArrayBox& flux,
            const amrex::Real* dx, double dt, amrex::Real gam, int Dim);


    std::vector<amrex_real> slopeLimiter(std::vector<double> delta_L, std::vector<double> delta_R,
            double omega, int Nv);

    void computeFluxForce(const amrex::FArrayBox& statein, const amrex::Box& bx,
        amrex::FArrayBox& flux, amrex::Real gam, const amrex::Real* dx, const amrex_real dt);

    void updateState(amrex::FArrayBox& stateout, const amrex::FArrayBox& statein,
             const amrex::Box& bx, const amrex::FArrayBox& flux, const amrex::Real* dx, const amrex_real dt, int Dim);

    void calcTimeStep(const amrex::FArrayBox& statein, const amrex::Real* dx, amrex::Real gam,
        amrex_real& dt);

    void convertConsArraytoPrim(amrex::FArrayBox& consIn,amrex::FArrayBox& primOut);
    void convertPrimArraytoCons(amrex::FArrayBox& primIn, amrex::FArrayBox& consOut);
    std::vector<amrex_real> Conserv2Prim(std::vector<amrex_real> U, amrex::Real gam);

    std::vector<amrex_real> Prim2Conserv(std::vector<amrex_real> W, amrex::Real gam);

    void linearAdv(const amrex::FArrayBox& statein, const amrex::Box& bx,
        amrex::FArrayBox& flux, amrex::Real gam, const amrex::Real* dx, const amrex::Real dt);

    void computeMaxGrad(const amrex::FArrayBox& statein, amrex::FArrayBox& rhoGradIn,
    const amrex::Box& bx, const amrex::Real* dx, amrex::Real& rhoGradMax);

    void performTaggingRhoGrad(amrex::Array4<char> tagfabArray, const amrex::Box& tilebx,
        const amrex::FArrayBox& statein,  const amrex::FArrayBox& rhoGradIn,
            amrex::Real rhoGradMax, const amrex::Real* dx);
    //Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;
    //
    //Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    //Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;

    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    //Do work after a restart().
    //
    virtual void post_restart () override;
    //
    //Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
			   int          n_error_buf = 0, int ngrow = 0) override;

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC.get(); }
#endif

    static int  NUM_STATE;
    static int  NUM_GROW;

protected:

    static void read_params ();

    // inline functioned implemented below
    MultiphaseAMR& getLevel (int lev);
    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The data.
    //
    amrex::FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static int          verbose;
    static amrex::Real  cfl;
    static int          do_reflux;
    static int          num_state_type;

#ifdef AMREX_PARTICLES
    void init_particles ();
    static int       do_tracers;
    static std::unique_ptr<amrex::AmrTracerParticleContainer> TracerPC;
#endif
};

//
// Inlines.
//

inline
MultiphaseAMR&
MultiphaseAMR::getLevel (int lev)
{
    return *(MultiphaseAMR *) &parent->getLevel(lev);
}

inline
amrex::FluxRegister&
MultiphaseAMR::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
MultiphaseAMR::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_MultiphaseAMR_H_*/
