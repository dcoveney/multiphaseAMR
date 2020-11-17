
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int      MultiphaseAMR::verbose         = 0;
Real     MultiphaseAMR::cfl             = 0.9;
int      MultiphaseAMR::do_reflux       = 0;

int      MultiphaseAMR::NUM_STATE       = 8;  // Four variables in the state
int      MultiphaseAMR::NUM_GROW        = 4;  // number of ghost cells
int      MultiphaseAMR::num_state_type  = 2;  // number of ghost cells

Real   gamL                          = 1.4;
Real   gamR                          = 1.4;
std::string probType                = "";
int    Ncons                           = 4;
int    Nprim                           = 4;
int    DEN                             = 0;
int    Xmom                            = 1;
int    Ymom                            = 2;
int    ENEtot                          = 3;
int    Xvel                            = 1;
int    Yvel                            = 2;
int    Press                           = 3;
int    ENEint                         = 4;

#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> MultiphaseAMR::TracerPC =  nullptr;
int MultiphaseAMR::do_tracers                       =  0;
#endif

//
//Default constructor.  Builds invalid object.
//
MultiphaseAMR::MultiphaseAMR ()
{
    flux_reg = 0;
}

//
//The basic constructor.
//
MultiphaseAMR::MultiphaseAMR (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

//
//The destructor.
//
MultiphaseAMR::~MultiphaseAMR ()
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
void
MultiphaseAMR::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void
MultiphaseAMR::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers and level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

//
//Write a plotfile to specified directory.
//
void
MultiphaseAMR::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);

#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}

//
//Define data descriptors.
//
void
MultiphaseAMR::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    // read_params();
    ParmParse pp_init("init");
    pp_init.query("probType", probType);
    desc_lst.addDescriptor(Prim_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,Ncons,
			   &cell_cons_interp);
    desc_lst.addDescriptor(Cons_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,Nprim,
			   &cell_cons_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	lo_bc[i] = hi_bc[i] = BCType::foextrap;   // transmissive boundaries
    }

    BCRec bc(lo_bc, hi_bc);

        desc_lst.setComponent(Cons_Type, DEN, "den", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Cons_Type, Xmom, "mom_x", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Cons_Type, Ymom, "mom_y", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Cons_Type, ENEtot, "E_t", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Prim_Type, DEN, "rho", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Prim_Type, Xvel, "u", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Prim_Type, Yvel, "v", bc,
                  StateDescriptor::BndryFunc(phifill));
        desc_lst.setComponent(Prim_Type, Press, "p", bc,
                  StateDescriptor::BndryFunc(phifill));

    // else if(probType == "ALLAIRE")
    // {
    //     desc_lst.setComponent(Prim, 0, "rho1", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 1, "rho2", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 2, "xmom", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 3, "ymom", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 4, "Et", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 5, "z", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 6, "uStar", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 7, "u", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 8, "v", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 9, "P", bc,
    //               StateDescriptor::BndryFunc(phifill));
    //     desc_lst.setComponent(Phi_Type, 10, "e", bc,
    //               StateDescriptor::BndryFunc(phifill));
    // }

}

//
//Cleanup data descriptors at end of run.
//
void
MultiphaseAMR::variableCleanUp ()
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif
}

//
//Initialize grid data at problem start-up.
//
void
MultiphaseAMR::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int testNumber;
    ParmParse pp("init");
    pp.query("testNumber", testNumber);
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Prim_Type);
    MultiFab& Cons_new = get_new_data(Cons_Type);
    Real cur_time   = state[Prim_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        RealBox    gridloc = RealBox(box,geom.CellSize(),geom.ProbLo());
        FArrayBox& Sfab = S_new[mfi];
        FArrayBox& Consfab = Cons_new[mfi];

        initialState(Sfab, box, gridloc, dx, testNumber);
        convertPrimArraytoCons(Sfab, Consfab);

          // initdata(&level, &cur_time, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
		  //  BL_TO_FORTRAN_3D(S_new[mfi]), &NUM_STATE, AMREX_ZFILL(dx),
		  //  AMREX_ZFILL(prob_lo));
    }

#ifdef AMREX_PARTICLES
    init_particles();
#endif

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level
                       << " data " << std::endl;
    }
}

//
//Initialize data on this level from another MultiphaseAMR (during regrid).
//
void
MultiphaseAMR::init (AmrLevel &old)
{
    MultiphaseAMR* oldlev = (MultiphaseAMR*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Cons_Type].curTime();
    Real prev_time = oldlev->state[Cons_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Cons_Type);

    FillPatch(old, S_new, 0, cur_time, Cons_Type, 0, Ncons);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
MultiphaseAMR::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Cons_Type].curTime();
    Real prev_time = getLevel(level-1).state[Cons_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    // MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_new = get_new_data(Cons_Type);
    FillCoarsePatch(S_new, 0, cur_time, Cons_Type, 0, Ncons);
}

//
//Advance grids at this level in time.
//
Real
MultiphaseAMR::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
    MultiFab& S_mm = get_new_data(Cons_Type);
    Real maxval = S_mm.max(0);
    Real minval = S_mm.min(0);
    std::ofstream amrexEulerOut;
    for (int k = 0; k < num_state_type; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    // MultiFab& S_new = get_new_data(Phi_Type);
    MultiFab& S_new = get_new_data(Cons_Type);
    MultiFab& Prim_new = get_new_data(Prim_Type);

    const Real prev_time = state[Cons_Type].prevTime();
    const Real cur_time = state[Cons_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
	fine = &getFluxReg(level+1);
	fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
	current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }


    // State with ghost cells
    MultiFab SBorder(grids, dmap, Ncons, NUM_GROW);
    MultiFab SReconstL(grids, dmap, Ncons, 3);
    MultiFab SReconstR(grids, dmap, Ncons, 3);
    FillPatch(*this, SBorder, NUM_GROW, time, Cons_Type, 0, Ncons);

    std::vector<double> UL(Ncons), UR(Ncons), UU(Ncons), UD(Ncons);
    std::vector<double> Fx(Ncons), Fy(Ncons);
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = SBorder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& ULbarx      =   SReconstL[mfi];
	    FArrayBox& URbarx      =   SReconstR[mfi];

	    // Allocate fabs for fluxes and Godunov velocities.
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bxtmp = amrex::surroundingNodes(bx,i);
		flux[i].resize(bxtmp,Ncons);
		uface[i].resize(amrex::grow(bxtmp, iteration), 1);
	    }
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        Array4<Real const> const& stateArray = statein.array();
        Array4<Real> const& fluxArray = flux[0].array();

        for(int i = lo.x; i <= hi.x+1; i++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int n = 0; n < Ncons; n++)
                    {
                        UL[n] = stateArray(i-1,j,k,n);
                        UR[n] = stateArray(i,j,k,n);
                    }

                    Fx = fluxHLLCcell(UR, UL, gam, 1, Ncons);
                    for(int n = 0; n < Ncons; n++)
                    {
                        fluxArray(i, j, k, n) = Fx[n];
                    }
                }
            }
        }
           // primMUSCLThirdOrder(statein, bx, ULbarx, URbarx, flux[0], dx, dt, gam, 1);
           updateState(stateout, statein, bx, flux[0], dx, dt, 1);

	    if (do_reflux) {
		    fluxes[0][mfi].copy(flux[0],mfi.nodaltilebox(0));
	    }
    }
    }
    if (do_reflux) {
    if (current) {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
    }
    // MultiFab SBorder(grids, dmap, NUM_STATE, NUM_GROW);

    MultiFab::Copy(SBorder, S_new, 0, 0, Ncons, S_new.nGrow());
    FillPatch(*this, SBorder, NUM_GROW, cur_time, Cons_Type, 0, Ncons);
    // FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = SBorder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& primout      =   Prim_new[mfi];
	    FArrayBox& ULbar      =   SReconstL[mfi];
	    FArrayBox& URbar      =   SReconstR[mfi];

	    // Allocate fabs for fluxes and Godunov velocities.
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bxtmp = amrex::surroundingNodes(bx,i);
		flux[i].resize(bxtmp,Ncons);
		uface[i].resize(amrex::grow(bxtmp, iteration), 1);
	    }
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        Array4<Real const> const& stateArray = statein.array();
        Array4<Real> const& fluxArray = flux[1].array();

        for(int i = lo.x; i <= hi.x; i++)
        {
            for(int j = lo.y; j <= hi.y+1; j++)
            {
                for(int k = lo.z; k <= hi.z; k++)
                {
                    for(int n = 0; n < Ncons; n++)
                    {
                        UD[n] = stateArray(i,j-1,k,n);
                        UU[n] = stateArray(i,j,k,n);
                    }

                    Fy = fluxHLLCcell(UU, UD, gam, 2, Ncons);
                    for(int n = 0; n < Ncons; n++)
                    {
                        fluxArray(i, j, k, n) = Fy[n];
                    }
                }
            }
        }
       // MUSCLHancock(statein, bx, ULbar, URbar, flux[1], dx, dt, gam, 2);
       updateState(stateout, statein, bx, flux[1], dx, dt, 2);

       convertConsArraytoPrim(stateout,primout);

	    if (do_reflux) {
		    fluxes[1][mfi].copy(flux[1],mfi.nodaltilebox(1));
	    }

    }
}
    if (do_reflux) {
	if (current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
	}
	if (fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
	}
    }


// #ifdef AMREX_PARTICLES
//     if (TracerPC) {
//       TracerPC->AdvectWithUmac(Umac, level, dt);
//     }
// #endif
    // amrex::Abort();
    return dt;

}

void
MultiphaseAMR::initialState(FArrayBox& statein, const Box& bx,
    RealBox gridloc, const Real* dx, int testNumber)
{
    Array4<Real> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    const Real* xlo = gridloc.lo();
    const Real* xhi = gridloc.hi();
    Real x, y, z;
    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                x = xlo[0] + dx[0]*((i-lo.x) + 0.5);
                y = xlo[1] + dx[1]*((j-lo.y) + 0.5);
                #if(AMREX_SPACEDIM==3)
                z = xlo[2] + dx[2]*((k-lo.z) + 0.5);
                #endif

                if(testNumber == 1)
                {
                    // Toro test 1
                    if(x < 0.5)
                    {
                        stateArray(i,j,k,DEN) = 1.0;
                        stateArray(i,j,k,Xvel) = 0.0;
                        stateArray(i,j,k,Yvel) = 0.0;
                        stateArray(i,j,k,Press) = 1.0;
                    }
                    else
                    {
                        stateArray(i,j,k,DEN) = 0.125;
                        stateArray(i,j,k,Xvel) = 0.0;
                        stateArray(i,j,k,Yvel) = 0.0;
                        stateArray(i,j,k,Press) = 0.1;
                    }
                }
            }
        }
    }
}
std::vector<amrex_real>
MultiphaseAMR::eqnFluxes(std::vector<amrex_real> U, amrex_real gam, int Dim, int Nv)
{
    std::vector<Real> F(Ncons);

	double p = (gam-1.0)*(U[ENEtot] - 0.5*(pow(U[Xmom],2.0)+pow(U[Ymom],2.0))/U[DEN]);

	if(Dim == 1)
	{
		// Equation fluxes in the x-direction.
		F[DEN] = U[Xmom];
		F[Xmom] = pow(U[Xmom],2.0)/U[DEN] + p;
		F[Ymom] = U[Xmom]*(U[Ymom]/U[DEN]);
		F[ENEtot] = (U[Xmom]/U[DEN])*(U[ENEtot]+p);
	}

	else if(Dim == 2)
	{
		// Equation fluxes in the y-direction.
		F[DEN] = U[Ymom];
		F[Xmom] = U[Ymom]*(U[Xmom]/U[DEN]);
		F[Ymom] = pow(U[Ymom],2.0)/U[DEN] + p;
		F[ENEtot] = (U[Ymom]/U[DEN])*(U[ENEtot]+p);
	}

	return F;
}

void
MultiphaseAMR::computeFluxForce(const FArrayBox& statein, const Box& bx,
    FArrayBox& flux, amrex_real gam, const amrex_real* dx, const amrex_real dt)
{

    Array4<Real const> const& stateArray = statein.array();
    Array4<Real> const& fluxArray_x = flux.array();
    // Array4<Real> const& fluxArray_y = flux.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    int Nv = stateArray.nComp();
    std::vector<Real> U_L(Nv);
    std::vector<Real> U_R(Nv);
    std::vector<Real> F_LF(Nv);
    std::vector<Real> U_mid(Nv);
    std::vector<Real> F_force(Nv);
    std::vector<Real> F_rm(Nv);
    std::vector<Real> F_l(Nv);
    std::vector<Real> F_r(Nv);
    Real c = 0.9;


    for(int k = lo.z; k <= hi.z; ++k)
    {
        for(int j = lo.y; j <= hi.y; ++j)
        {
            for(int i = lo.x; i <= hi.x+1; ++i)
            {
                for(int n = 0; n < Nv; n++)
                {
                        U_L[n] = stateArray(i-1,j,k,n);
                        U_R[n] = stateArray(i,j,k,n);
                }


                F_l = eqnFluxes(U_L, gam, 1, Nv);
                F_r = eqnFluxes(U_R, gam, 1, Nv);
                for ( int p = 0 ; p < Nv ; p++ )
                {
                     //F_LF[p] = ((1+c)/(2*c))*F_l[p] + ((c-1)/(2*c))*F_r[p];
                    F_LF[p] = 0.5*(F_l[p] + F_r[p]) - (dx[0]/dt)*0.5*(U_R[p] - U_L[p]);

                    U_mid[p] = 0.5*(U_L[p] + U_R[p]) - 0.5*(dt/dx[0])*(F_r[p] - F_l[p]);
                }

                F_rm = eqnFluxes(U_mid, gam, 1, Nv);
                //
                // for(int p = 0; p < Nv; p++)
                // {
                //     F_force[p] = 0.5*(F_LF[p]+F_rm[p]);
                // }

                F_force = fluxHLLCcell(U_R, U_L, gam, 1, Nv);
                for(int n = 0; n < Nv; n++)
                {
                    fluxArray_x(i,j,k,n) = F_force[n];
                }
                // if(j == 5 || k == 5)
                // {
                //     std::cout << fluxArray_x(i, j, k, 0) << " ";
                // }
            }
        }

    }
//     std::cout << "Fluxes," << std::endl;
//
//     for(int i = lo.x; i <= hi.x+1; i++)
//     {
//             std::cout << fluxArray_x(i,0,0,0) << " ";
//
//     }
//
//
// std::cout << std::endl;

    //
    // F_l = eqnFluxes(U_L, gam, 1);
    // F_r = eqnFluxes(U_R, gam, 1);
    //
    // for ( int i = 0 ; i < Nv ; i++ )
	// {
	// 	F_LF[i] = 0.5*(F_l[i] + F_r[i]) - (dx[0]/dt)*0.5*(U_R[i] - U_L[i]);
    //
	// 	U_mid[i] = 0.5*(U_L[i] + U_R[i]) - 0.5*(dt/dx[0])*(F_r[i] - F_l[i]);
	// }
    //
    // F_rm = eqnFluxes(U_mid, gam, 1);
    //
    // for(int i = 0; i < Nv; i++)
    // {
    //     F_force[i] = 0.5*(F_LF[i] + F_rm[i]);
    // }
    //
    // for(int k = lo.z; k <= hi.z; k++)
    // {
    //     for(int j = lo.y; j <= hi.y; j++)
    //     {
    //         for(int i = lo.x; i <= hi.x; i++)
    //         {
    //             for(int n = 0; n < Nv; n++)
    //             {
    //                 fluxArray_x(i,j,k,n) = F_force[n];
    //             }
    //         }
    //     }
    // }
    //

}

std::vector<Real>
MultiphaseAMR::slopeLimiter(std::vector<double> delta_L, std::vector<double> delta_R, double omega,
    int Nv)
{

	double r_R, r_L;
	double phi, phi_R, phi_L;
	std::vector<double> r(Nv);
	std::vector<double> xi_L(Nv);
	std::vector<double> xi_R(Nv);
	std::vector<double> xi(Nv);

	double beta_L = 1; // 2.0/(1+CFL);
	double beta_R = 1; //2.0/(1-CFL);

	double minDummy;
	// 0 for superbee, 1 for minbee, 2 for van-leer
	int limiter_fcn = 2;

	for(int i = 0; i < Nv; i++)
	{
		r[i] = delta_L[i]/delta_R[i];
		xi_L[i] = (2.0*beta_L*r[i])/(1.0-omega + (1.0+omega)*r[i]);
		xi_R[i] = (2.0*beta_R)/(1.0-omega + (1.0+omega)*r[i]);


		if(delta_L[i] < 0.0 && delta_R[i] > 0.0)
		{
			xi[i] = 0.0;
		}
		else if(delta_L[i] == 0.0 && delta_R[i] == 0.0)
		{
			xi[i] = 1.0;
		}

		else if(delta_L[i] != 0.0 && delta_R[i] == 0.0)
		{
			xi[i] = 0.0;
		}

		else if(isnan(r[i]))
		{
			xi[i] = 1.0;
		}
		else if(isinf(r[i]))
		{
			xi[i] = 0.0;
		}

		if(limiter_fcn == 0)
		{
			if(r[i] < 0.0)
			{
				xi[i] = 0.0;
			}

			if(r[i] >= 0.0 && r[i] < 0.5)
			{
				xi[i] = 2.0*r[i];
			}
			if(r[i] >= 0.5 && r[i] < 1.0)
			{
				xi[i] = 1.0;
			}
			if(r[i] >= 1.0)
			{
				minDummy = min(r[i], xi_R[i]);
				xi[i] = min(minDummy, 2.0);
			}
		}

		if(limiter_fcn == 1)
		{
			if(r[i] < 0.0)
			{
				xi[i] = 0.0;
			}
			if(r[i] >= 0.0 && r[i] < 1.0)
			{
				xi[i] = r[i];
			}
			if(r[i] >= 1.0)
			{
				xi[i] = min(1.0, xi_R[i]);
			}
		}

		if(limiter_fcn == 2)
		{
			if(r[i] <= 0.0)
			{
				xi[i] = 0.0;
			}
			if(r[i] > 0.0)
			{
				xi[i] = fmin((2.0*r[i])/(1.0+r[i]), xi_R[i]);
			}
		}
	}

	return xi;

}
//
//
void
MultiphaseAMR::MUSCLHancock(const FArrayBox& statein, const Box& bx,
    FArrayBox& ULbar, FArrayBox& URbar, FArrayBox& flux,
        const amrex_real* dx, double dt, amrex_real gam, int Dim)
{
    Array4<Real const> const& stateArray = statein.array();
    Array4<Real> const& fluxArray = flux.array();
    int Nv = 4;
	std::vector<double> U_L(Nv);
	std::vector<double> U_LL(Nv);
	std::vector<double> U_R(Nv);
	std::vector<double> U_RR(Nv);
	std::vector<double> Ul(Nv);
	std::vector<double> Ur(Nv);
	std::vector<double> F_l(Nv);
	std::vector<double> F_r(Nv);
	// vector<vector<double>> delta_Rvec;
	// vector<vector<double>> delta_Lvec;
	std::vector<double> delta_L(Nv);
	std::vector<double> delta_R(Nv);
	std::vector<double> delta(Nv);
	std::vector<double> deltaLim(Nv);
	std::vector<double> xi(Nv);
    // Array4<Real> const& fluxArray_y = flux.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    std::vector<double> U_Rbar_cell(Nv);
    std::vector<double> U_Lbar_cell(Nv);
    std::vector<double> F_hllc;
    // std::vector<std::vector<double>> U_Lbar;
    // std::vector<std::vector<double>> U_Rbar;
    // U_Lbar.resize(Nv, std::vector<double>(hi.x+2));
    // U_Rbar.resize(Nv, std::vector<double>(hi.x+2));
    Array4<Real> const& U_Lbar = ULbar.array();
    Array4<Real> const& U_Rbar = URbar.array();
	double omega = 0.0;

    if(Dim == 1)
    {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif

        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int i = lo.x; i <= hi.x+2; i++)
                {
                    for(int n = 0; n < Nv; n++)
            		{
            			 U_LL[n] = stateArray(i-2, j, k, n);
            			 U_L[n] =  stateArray(i-1, j, k, n);
            			 U_R[n] =  stateArray(i, j, k, n);
            	 	}

            		 for(int n = 0; n < Nv; n++)
            		 {
            			 // delta_Lvec[j][i] = U_L[j] - U_LL[j];
            			 // delta_Rvec[j][i] = U_R[j] - U_L[j];
            			 delta_L[n] = U_L[n] - U_LL[n];
            			 delta_R[n] = U_R[n] - U_L[n];
            		 }

            		 xi = slopeLimiter(delta_L, delta_R, omega, Nv);

            		 for(int n = 0; n < Nv; n++)
            		 {
            			 delta[n] = 0.5*(1.0+omega)*delta_L[n] + 0.5*(1.0-omega)*delta_R[n];
            			 deltaLim[n] =  xi[n]*delta[n];
            			 Ul[n] = U_L[n] - 0.5*deltaLim[n];
            			 Ur[n] = U_L[n] + 0.5*deltaLim[n];
            		 }

            		 F_l = eqnFluxes(Ul, gam, Dim, Nv);
            	 	 F_r = eqnFluxes(Ur, gam, Dim, Nv);

             		 for(int n = 0; n < Nv; n++)
             		 {
             			 U_Lbar(i,j,k,n) = Ul[n] + 0.5*(dt/dx[0])*(F_l[n] - F_r[n]);
             			 U_Rbar(i,j,k,n) = Ur[n] + 0.5*(dt/dx[0])*(F_l[n] - F_r[n]);
             		 }
            	 }
             }
         }

         #ifdef _OPENMP
         #pragma omp parallel for
         #endif
         for(int k = lo.z; k <= hi.z; k++)
         {
             for(int j = lo.y; j <= hi.y; j++)
             {
                 for(int i = lo.x; i <= hi.x+1; i++)
                 {
                     for(int n = 0; n < Nv; n++)
             		{
                        U_Rbar_cell[n] = U_Rbar(i,j,k,n);
                        U_Lbar_cell[n] = U_Lbar(i+1,j,k,n);
                    }

                    F_hllc = fluxHLLCcell(U_Lbar_cell, U_Rbar_cell, gam, Dim, Nv);
                    for(int n = 0; n < Nv; n++)
                    {
                        fluxArray(i, j, k, n) = F_hllc[n];
                    }
                }
            }
        }

    }

    else if(Dim == 2)
    {
        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y+2; j++)
            {
                for(int i = lo.x; i <= hi.x; i++)
                {
                    for(int n = 0; n < Nv; n++)
                    {
                         U_LL[n] = stateArray(i, j-2, k, n);
                         U_L[n] =  stateArray(i, j-1, k, n);
                         U_R[n] =  stateArray(i, j, k, n);
                    }

                     for(int n = 0; n < Nv; n++)
                     {
                         // delta_Lvec[j][i] = U_L[j] - U_LL[j];
                         // delta_Rvec[j][i] = U_R[j] - U_L[j];
                         delta_L[n] = U_L[n] - U_LL[n];
                         delta_R[n] = U_R[n] - U_L[n];
                     }

                     xi = slopeLimiter(delta_L, delta_R, omega, Nv);

                     for(int n = 0; n < Nv; n++)
                     {
                         delta[n] = 0.5*(1.0+omega)*delta_L[n] + 0.5*(1.0-omega)*delta_R[n];
                         deltaLim[n] = xi[n]*delta[n];
                         Ul[n] = U_L[n] - 0.5*deltaLim[n];
                         Ur[n] = U_L[n] + 0.5*deltaLim[n];
                     }

                     F_l = eqnFluxes(Ul, gam, Dim, Nv);
                     F_r = eqnFluxes(Ur, gam, Dim, Nv);

                     for(int n = 0; n < Nv; n++)
                     {
                         U_Lbar(i,j,k,n) = Ul[n] + 0.5*(dt/dx[1])*(F_l[n] - F_r[n]);
                         U_Rbar(i,j,k,n) = Ur[n] + 0.5*(dt/dx[1])*(F_l[n] - F_r[n]);
                     }
                 }
             }
         }

         for(int k = lo.z; k <= hi.z; k++)
         {
             for(int j = lo.y; j <= hi.y+1; j++)
             {
                 for(int i = lo.x; i <= hi.x; i++)
                 {
                     for(int n = 0; n < Nv; n++)
                    {
                        U_Rbar_cell[n] = U_Rbar(i,j,k,n);
                        U_Lbar_cell[n] = U_Lbar(i,j+1,k,n);
                    }

                    F_hllc = fluxHLLCcell(U_Lbar_cell, U_Rbar_cell, gam, Dim, Nv);
                    for(int n = 0; n < Nv; n++)
                    {
                        fluxArray(i, j, k, n) = F_hllc[n];
                    }
                }
            }
        }
    }


 }

 void MultiphaseAMR::primMUSCLThirdOrder(const FArrayBox& primIn, const Box& bx,
     FArrayBox& ULbar, FArrayBox& URbar, FArrayBox& flux,
         const amrex_real* dx, double dt, amrex_real gam, int Dim)
 {
     Array4<Real const> const& primArr = primIn.array();
     Array4<Real> const& fluxArray = flux.array();

     std::vector<double> W_L(Nprim);
     std::vector<double> W_LL(Nprim);
     std::vector<double> W_R(Nprim);
     std::vector<double> W_RR(Nprim);
     std::vector<double> Wl(Nprim);
     std::vector<double> Wr(Nprim);
     std::vector<double> Ul(Nprim);
     std::vector<double> Ur(Nprim);
     std::vector<double> F_l(Nprim);
     std::vector<double> F_r(Nprim);
     std::vector<double> F_hllc(Nprim);
     Dim3 lo = lbound(bx);
     Dim3 hi = ubound(bx);
     double kap = 1.0/3.0;
     std::vector<double> delta_L(Nprim);
     std::vector<double> delta_R(Nprim);
     std::vector<double> delta(Nprim);
     std::vector<double> rR(Nprim), rL(Nprim), phiL(Nprim), phiR(Nprim);
     Array4<Real> const& U_Lbar = ULbar.array();
     Array4<Real> const& U_Rbar = URbar.array();
     std::vector<double> U_Rbar_cell(Nprim);
     std::vector<double> U_Lbar_cell(Nprim);
     if(Dim == 1)
     {
         #ifdef _OPENMP
         #pragma omp parallel for
         #endif

         for(int k = lo.z; k <= hi.z; k++)
         {
             for(int j = lo.y; j <= hi.y; j++)
             {
                 for(int i = lo.x; i <= hi.x+1; i++)
                 {
                     for(int n = 0; n < Nprim; n++)
             		{
             			 W_LL[n] = primArr(i-2, j, k, n);
             			 W_L[n] =  primArr(i-1, j, k, n);
             			 W_R[n] =  primArr(i, j, k, n);
             			 W_RR[n] =  primArr(i+1, j, k, n);
             			 // delta_Lvec[j][i] = U_L[j] - U_LL[j];
             			 // delta_Rvec[j][i] = U_R[j] - U_L[j];
             			 delta_L[n] = W_L[n] - W_LL[n];
             			 delta[n] = W_R[n] - W_L[n];
             			 delta_R[n] = W_RR[n] - W_R[n];
                         rR[n] = delta[n]/delta_R[n];
                         rL[n] = delta_L[n]/delta[n];


                         phiL[n] = fmax(0.0,(2.0*rL[n])/(rL[n]*rL[n]+1.0));
                         phiR[n] = fmax(0.0,(2.0*rR[n])/(rR[n]*rR[n]+1.0));

                         if(isnan(rR[n]))
                       {
                           phiR[n] = 1.0;
                       }
                       else if(isinf(rR[n]))
                       {
                           phiR[n] = 0.0;
                       }
                         if(isnan(rL[n]))
                       {
                           phiL[n] = 1.0;
                       }
                       else if(isinf(rL[n]))
                       {
                           phiL[n] = 0.0;
                       }

                         Wl[n] = W_L[n] + (phiL[n]/4.0)*((1.0-kap*phiL[n])*delta_L[n]+(1.0+kap*phiL[n])*delta[n]);
                         Wr[n] = W_R[n] + (phiR[n]/4.0)*((1.0-kap*phiR[n])*delta_R[n]+(1.0+kap*phiR[n])*delta[n]);
             		 }

                     Ul = Prim2Conserv(Wl,gam);
                     Ur = Prim2Conserv(Wr,gam);
                     F_l = eqnFluxes(Ul, gam, Dim, Nprim);
                     F_r = eqnFluxes(Ur, gam, Dim, Nprim);

                     for(int n = 0; n < Nprim; n++)
                     {
                         U_Lbar(i,j,k,n) = Ul[n] + 0.5*(dt/dx[0])*(F_l[n] - F_r[n]);
                         U_Rbar(i,j,k,n) = Ur[n] + 0.5*(dt/dx[0])*(F_l[n] - F_r[n]);
                     }
                 }
             }
         }
         #ifdef _OPENMP
         #pragma omp parallel for
         #endif
         for(int k = lo.z; k <= hi.z; k++)
         {
             for(int j = lo.y; j <= hi.y; j++)
             {
                 for(int i = lo.x; i <= hi.x+1; i++)
                 {
                     for(int n = 0; n < Nprim; n++)
                   {
                        U_Rbar_cell[n] = U_Rbar(i,j,k,n);
                        U_Lbar_cell[n] = U_Lbar(i+1,j,k,n);
                    }

                    F_hllc = fluxHLLCcell(U_Lbar_cell, U_Rbar_cell, gam, Dim, Ncons);
                    for(int n = 0; n < Nprim; n++)
                    {
                        fluxArray(i, j, k, n) = F_hllc[n];
                    }
                }
            }
        }
     }


 }
//
std::vector<amrex_real>
MultiphaseAMR::fluxHLLCcell(std::vector<double>& U_R, std::vector<double>& U_L,
    amrex_real gam, int Dim, int Nv)
{
    bool speed_estimate = 1;

	std::vector<Real> W_L(Nv);
	std::vector<Real> W_R(Nv);

	double S_L, S_R;

	std::vector<Real> F_L(Nv);
	std::vector<Real> F_R(Nv);

	std::vector<Real> F_Lstar(Nv);
	std::vector<Real> F_Rstar(Nv);

	std::vector<Real> U_Lstar(Nv);
	std::vector<Real> U_Rstar(Nv);

	double S_Star;

	double K_L, K_R;

	std::vector<Real> F_hllc(Nv);


	 W_L = Conserv2Prim(U_L, gam);
	 W_R = Conserv2Prim(U_R, gam);

	 if(Dim == 1)
	 {

		 if(speed_estimate == 0)
		 {

			 // Speed estimate based on Roe averages:

			 Real H_L, H_R;
			 Real u_tilde, H_tilde, c_tilde;

			 H_L = (U_L[3] + W_L[3])/U_L[0];
			 H_R = (U_R[3] + W_R[3])/U_R[0];

			 u_tilde = (sqrt(U_L[0])*W_L[1] + sqrt(U_R[0])*W_R[1])/(sqrt(U_L[0]) + sqrt(U_R[0]));

			 H_tilde = (sqrt(U_L[0])*H_L + sqrt(U_R[0])*H_R)/(sqrt(U_L[0]) + sqrt(U_R[0]));

			 c_tilde = sqrt((gam-1)*(H_tilde - 0.5*pow(u_tilde, 2)));

			 S_L = u_tilde - c_tilde;
			 S_R = u_tilde + c_tilde;
		 }

		 else
		 {
			 S_L = fmin(W_L[1] - sqrt(gam*W_L[3]/W_L[0]), W_R[1] - sqrt(gam*W_R[3]/W_R[0]));
			 S_R = fmax(W_L[1] + sqrt(gam*W_L[3]/W_L[0]), W_R[1] + sqrt(gam*W_R[3]/W_R[0]));
		 }



		 S_Star = (W_R[3] - W_L[3] + U_L[1]*(S_L - W_L[1]) - U_R[1]*(S_R - W_R[1])) /
		 (W_L[0]*(S_L-W_L[1]) - W_R[0]*(S_R - W_R[1]));

		 K_L = W_L[0]*(S_L - W_L[1])/(S_L - S_Star);
		 K_R = W_R[0]*(S_R - W_R[1])/(S_R - S_Star);

		 U_Lstar[0] = K_L;
		 U_Lstar[1] = K_L*S_Star;
		 U_Lstar[2] = K_L*(W_L[2]);
		 U_Lstar[3] = K_L*(U_L[3]/U_L[0] + (S_Star - W_L[1])*(S_Star + W_L[3]/(W_L[0]*(S_L - W_L[1]))));

		 U_Rstar[0] = K_R;
		 U_Rstar[1] = K_R*S_Star;
		 U_Rstar[2] = K_R*(W_R[2]);
		 U_Rstar[3] = K_R*(U_R[3]/U_R[0] + (S_Star - W_R[1])*(S_Star + W_R[3]/(W_R[0]*(S_R - W_R[1]))));

		 F_L = eqnFluxes(U_L, gam, Dim, Nv);
		 F_R = eqnFluxes(U_R, gam, Dim, Nv);

		 for(int i = 0; i < Nv; i++)
		 {
			 F_Lstar[i] = F_L[i] + S_L*(U_Lstar[i] - U_L[i]);
			 F_Rstar[i] = F_R[i] + S_R*(U_Rstar[i] - U_R[i]);

		 }

		 if(S_L >= 0)
		 {
			 F_hllc = F_L;
		 }

		 if( S_L < 0 && S_Star >= 0)
		 {
			 F_hllc = F_Lstar;
		 }

		 if(S_Star < 0 && S_R >= 0)
		 {
			 F_hllc = F_Rstar;
		 }

		 if(S_R < 0)
		 {
			 F_hllc = F_R;
		 }

	 }

	 else if(Dim == 2)
	 {
		 if(speed_estimate == 0)
		 {

			 // Speed estimate based on Roe averages:

			 Real H_L, H_R;
			 Real u_tilde, H_tilde, c_tilde;

			 H_L = (U_L[3] + W_L[3])/U_L[0];
			 H_R = (U_R[3] + W_R[3])/U_R[0];

			 u_tilde = (sqrt(U_L[0])*W_L[2] + sqrt(U_R[0])*W_R[2])/(sqrt(U_L[0]) + sqrt(U_R[0]));

			 H_tilde = (sqrt(U_L[0])*H_L + sqrt(U_R[0])*H_R)/(sqrt(U_L[0]) + sqrt(U_R[0]));

			 c_tilde = sqrt((gam-1)*(H_tilde - 0.5*pow(u_tilde, 2)));

			 S_L = u_tilde - c_tilde;
			 S_R = u_tilde + c_tilde;
		 }

		 else
		 {
			 S_L = fmin(W_L[2] - sqrt(gam*W_L[3]/W_L[0]), W_R[2] - sqrt(gam*W_R[3]/W_R[0]));
			 S_R = fmax(W_L[2] + sqrt(gam*W_L[3]/W_L[0]), W_R[2] + sqrt(gam*W_R[3]/W_R[0]));
		 }



		 S_Star = (W_R[3] - W_L[3] + U_L[2]*(S_L - W_L[2]) - U_R[2]*(S_R - W_R[2])) /
		 (W_L[0]*(S_L-W_L[2]) - W_R[0]*(S_R - W_R[2]));

		 K_L = W_L[0]*(S_L - W_L[2])/(S_L - S_Star);
		 K_R = W_R[0]*(S_R - W_R[2])/(S_R - S_Star);

		 U_Lstar[0] = K_L;
		 U_Lstar[1] = K_L*(U_L[1]/U_L[0]);
		 U_Lstar[2] = K_L*S_Star;
		 U_Lstar[3] = K_L*(U_L[3]/U_L[0] + (S_Star - W_L[2])*(S_Star + W_L[3]/(W_L[0]*(S_L - W_L[2]))));

		 U_Rstar[0] = K_R;
		 U_Rstar[1] = K_R*(U_R[1]/U_R[0]);
		 U_Rstar[2] = K_R*S_Star;
		 U_Rstar[3] = K_R*(U_R[3]/U_R[0] + (S_Star - W_R[2])*(S_Star + W_R[3]/(W_R[0]*(S_R - W_R[2]))));


		 F_L = eqnFluxes(U_L, gam, Dim, Nv);
		 F_R = eqnFluxes(U_R, gam, Dim, Nv);

		 for(int i = 0; i < Nv; i++)
		 {
			 F_Lstar[i] = F_L[i] + S_L*(U_Lstar[i] - U_L[i]);
			 F_Rstar[i] = F_R[i] + S_R*(U_Rstar[i] - U_R[i]);

		 }

		 if(S_L >= 0)
		 {
			 F_hllc = F_L;
		 }

		 if( S_L < 0 && S_Star >= 0)
		 {
			 F_hllc = F_Lstar;
		 }

		 if(S_Star < 0 && S_R >= 0)
		 {
			 F_hllc = F_Rstar;
		 }

		 if(S_R < 0)
		 {
			 F_hllc = F_R;
		 }

	 }

	return F_hllc;
}


void
MultiphaseAMR::linearAdv(const FArrayBox& statein, const Box& bx,
    FArrayBox& flux, double gam, const amrex_real* dx, const amrex_real dt)
{
    Array4<Real const> const& stateArray = statein.array();
    Array4<Real> const& fluxArray_x = flux.array();
    // Array4<Real> const& fluxArray_y = flux.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    int Nv = stateArray.nComp();
    std::vector<Real> U_L(Nv);
    std::vector<Real> U_R(Nv);
    std::vector<Real> F(Nv);
    std::vector<Real> F_LF(Nv);
    std::vector<Real> U_mid(Nv);
    std::vector<Real> F_force(Nv);
    std::vector<Real> F_rm(Nv);
    Real advSpeed = 1.0;
    Real c = (advSpeed*dt)/dx[0];


    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x+1; i++)
            {
                for(int n = 0; n < Nv; n++)
                {
                        U_L[n] = stateArray(i-1,j,k,n);
                        U_R[n] = stateArray(i,j,k,n);

                        F_LF[n] = ((1+c)/(2*c))*(advSpeed*U_L[n]) + ((c-1)/(2*c))*(advSpeed*U_R[n]);

                        U_mid[n] = 0.5*(U_L[n] + U_R[n]) - 0.5*(dt/dx[0])*advSpeed*(U_R[n] - U_L[n]);
                        F_rm[n] = advSpeed*U_mid[n];

                        fluxArray_x(i,j,k,n) = 0.5*(F_LF[n]+F_rm[n]);
                }
            }
        }
    }


}

void
MultiphaseAMR::updateState(FArrayBox& stateout, const FArrayBox& statein,
     const Box& bx, const FArrayBox& flux, const amrex_real* dx, const amrex_real dt, int Dim)
{
    Array4<Real const> const& stateOld = statein.array();
    Array4<Real> const& stateNew = stateout.array();
    Array4<Real const> const& fluxArray = flux.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    int Nv = 4;


    if(Dim == 1)
    {
        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int i = lo.x; i <= hi.x; i++)
                {
                    for(int n = 0; n < Nv; n++)
                    {
                        // stateNew(i,j,k,n) = fluxArray(i,j,k,n)-fluxArray(i+1,j,k,n);
                        stateNew(i,j,k,n) = stateOld(i,j,k,n) +
                            (dt/dx[0])*(fluxArray(i,j,k,n)-fluxArray(i+1,j,k,n));
                    }
                    // std::cout << fluxArray_x(i+1,j,0,0) << " ";

                }
                // std::cout << std::endl;
            }
        }
    }

    else if(Dim == 2)
    {
        for(int k = lo.z; k <= hi.z; k++)
        {
            for(int j = lo.y; j <= hi.y; j++)
            {
                for(int i = lo.x; i <= hi.x; i++)
                {
                    for(int n = 0; n < Nv; n++)
                    {
                        // stateNew(i,j,k,n) = fluxArray(i,j,k,n)-fluxArray(i+1,j,k,n);
                        stateNew(i,j,k,n) = stateOld(i,j,k,n) +
                            (dt/dx[1])*(fluxArray(i,j,k,n)-fluxArray(i,j+1,k,n));
                    }
                    // std::cout << fluxArray_x(i+1,j,0,0) << " ";

                }
                // std::cout << std::endl;
            }
        }
    }

    // std::cout << "State out," << std::endl;
    //
    // for(int i = lo.x; i <= hi.x; i++)
    // {
    //         std::cout << stateNew(i,0,0,0) << " ";
    //
    // }


// std::cout << std::endl;

}

std::vector<amrex_real>
MultiphaseAMR::Conserv2Prim(std::vector<amrex_real> U, amrex_real gam)
{
    std::vector<Real> W(Nprim);
    W[DEN] = U[DEN];
    W[Xvel] = U[Xmom]/U[DEN];
    W[Yvel] = U[Ymom]/U[DEN];
    W[Press] = (gam - 1.0)*(U[ENEtot] - 0.5*(pow(U[Xmom],2.0)+pow(U[Ymom], 2.0))/U[DEN]);

    return W;
}

std::vector<amrex_real>
MultiphaseAMR::Prim2Conserv(std::vector<amrex_real> W, amrex_real gam)
{
	std::vector<double> U(4);
	U[DEN] = W[DEN];
	U[Xmom] = W[DEN]*W[Xvel];
	U[Ymom] = W[DEN]*W[Yvel];
	U[ENEtot] = W[Press]/(gam - 1.0) + 0.5*(W[DEN])*(pow(W[Xvel],2.0) + pow(W[Yvel], 2.0));

	return U;
}

void
MultiphaseAMR::calcTimeStep(const FArrayBox& statein, const Real* dx, amrex_real gam,
    amrex_real& dt)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(stateArray);
    Dim3 hi = ubound(stateArray);
    int Nv = stateArray.nComp();
    std::vector<Real> WCell(Nv);
    std::vector<Real> UCell(Nv);

    Real cS;
    Real SMax_x = 0.0;
    Real SMax_y = 0.0;
    Real CFL = 0.9;
    Real advSpeed = 1.0;
    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                for(int n = 0; n < Nprim; n++)
                {
                    UCell[n] = stateArray(i,j,k,n);
                }

                WCell = Conserv2Prim(UCell, gam);


                cS = sqrt(gam*WCell[Press]/WCell[DEN]);
                SMax_x = fmax(SMax_x, cS + fabs(WCell[Xvel]));
    			SMax_y = fmax(SMax_y, cS + fabs(WCell[Yvel]));
            }
        }
    }

     // dt = CFL*(dx[0]);
    dt = CFL*fmin(dx[0]/SMax_x, dx[1]/SMax_y);
    if(dt <= 1.0e-15)
    {
        amrex::Abort("Time step too small!");
    }
    // std::cout << dt << std::endl;
}

void
MultiphaseAMR::convertConsArraytoPrim(FArrayBox& consIn, FArrayBox& primOut)
{
    Array4<Real> const& primArr = primOut.array();
    Array4<Real> const& consArr = consIn.array();
    Dim3 lo = lbound(primArr);
    Dim3 hi = ubound(primArr);
    int Nv = primArr.nComp();

    std::vector<Real> consVars(Ncons);
    std::vector<Real> primVars(Nprim);

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {
                consVars[0] = consArr(i,j,k,DEN);
                consVars[1] = consArr(i,j,k,Xmom);
                consVars[2] = consArr(i,j,k,Ymom);
                consVars[3] = consArr(i,j,k,ENEtot);
                primVars = Conserv2Prim(consVars, gam);

                primArr(i,j,k,DEN) = primVars[0];
                primArr(i,j,k,Xvel) = primVars[1];
                primArr(i,j,k,Yvel) = primVars[2];
                primArr(i,j,k,Press) = primVars[3];
            }
        }
    }
}
void
MultiphaseAMR::convertPrimArraytoCons(FArrayBox& primIn, FArrayBox& consOut)
{
    Array4<Real> const& primArr = primIn.array();
    Array4<Real> const& consArr = consOut.array();
    Dim3 lo = lbound(primArr);
    Dim3 hi = ubound(primArr);
    int Nv = primArr.nComp();

    std::vector<Real> consVars(Nv);
    std::vector<Real> primVars(Nv);

    for(int k = lo.z; k <= hi.z; k++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int i = lo.x; i <= hi.x; i++)
            {

                primVars[0] = primArr(i,j,k,DEN);
                primVars[1] = primArr(i,j,k,Xvel);
                primVars[2] = primArr(i,j,k,Yvel);
                primVars[3] = primArr(i,j,k,Press);

                consVars = Prim2Conserv(primVars, gam);

                consArr(i,j,k,DEN) = consVars[0];
                consArr(i,j,k,Xmom) = consVars[1];
                consArr(i,j,k,Ymom) = consVars[2];
                consArr(i,j,k,ENEtot) = consVars[3];
                // stateArray(i,j,k,7) = primVars[3]/(primVars[0]*(gam-1));
            }
        }
    }
}

//Estimate time step.
//
Real
MultiphaseAMR::estTimeStep (Real)
{
    // This is just a dummy value to start with
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    // std::cout << dx[1] << std::endl;
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Cons_Type].curTime();
    MultiFab& S_new = get_new_data(Cons_Type);
    MultiFab SBorder(grids, dmap, Ncons, NUM_GROW);
    MultiFab::Copy(SBorder, S_new, 0, 0, Ncons, S_new.nGrow());
    FillPatch(*this, SBorder, NUM_GROW, cur_time, Cons_Type, 0, Ncons);
    Real dtCycle = 0.0;
#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	FArrayBox uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();
        //FArrayBox& statein = Cons_new[mfi];
        // Array4<Real> const& stateArray = statein.array();
        //     Dim3 lo = lbound(stateArray);
        //     Dim3 hi = ubound(stateArray);
        //     int Nv = statein.nComp();

            // for(int n = 0; n < Nv; n++)
            // {
            //     for(int k = lo.z; k <= hi.z; k++)
            //     {
            //         for(int j = lo.y; j <= hi.y; j++)
            //         {
            //             for(int i = lo.x; i <= hi.x; i++)
            //             {
            //                 std::cout << stateArray(i,j,k,n) << " ";
            //             }
            //         }
            //         std::cout << std::endl;
            //     }
            //
         //    // }
         // for (int i = 0; i < BL_SPACEDIM ; i++) {
         //     const Box& bx = mfi.nodaltilebox(i);
		 //     // uface[i].resize(bx,1);
         //     }
            const FArrayBox& statein = SBorder[mfi];
            calcTimeStep(statein, dx, gam, dtCycle);
            dt_est = fmin(dt_est, dtCycle);
        // computeTimeStep(BL_TO_FORTRAN_3D(statein),
        //     bx.loVect(), bx.hiVect(), &NUM_STATE, dx, &gam, &dt_est);
	    // for (int i = 0; i < BL_SPACEDIM; ++i) {
		// Real umax = uface[i].norm(0);
		// if (umax > 1.e-100) {
		//     dt_est = std::min(dt_est, dx[i] / umax);
		// }
	    // }
	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    // dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "MultiphaseAMR::estTimeStep at level " << level
                       << ":  dt_est = " << dt_est << std::endl;
    }

    return dt_est;
}

//
//Compute initial time step.
//
Real
MultiphaseAMR::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
MultiphaseAMR::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Cons_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
MultiphaseAMR::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        MultiphaseAMR& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1)
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Cons_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
MultiphaseAMR::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);

        if (iteration < ncycle || level == 0)
	  {
            int ngrow = (level == 0) ? 0 : iteration;

	    TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
	  }
      }
#endif
}

//
//Do work after regrid().
//
void
MultiphaseAMR::post_regrid (int lbase, int new_finest) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#endif
}

//
//Do work after a restart().
//
void
MultiphaseAMR::post_restart()
{
#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

//
//Do work after init().
//
void
MultiphaseAMR::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
MultiphaseAMR::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();
    Real rhoGradMaxCycle = 0.0;
    Real rhoGradMax = 0.0;
    const Real cur_time = state[Cons_Type].curTime();
    MultiFab& S_new = get_new_data(Cons_Type);
    MultiFab SBorder(grids, dmap, Ncons, NUM_GROW);
    MultiFab::Copy(SBorder, S_new, 0, 0, Ncons, S_new.nGrow());
    FillPatch(*this, SBorder, NUM_GROW, cur_time, Cons_Type, 0, Ncons);
    MultiFab rhoGradFab(grids, dmap, 1, NUM_GROW);


#ifdef _OPENMP
#pragma omp parallel
// #pragma omp parallel reduction(min:dt_est)
#endif
    {
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
        const Box& bx = mfi.tilebox();
        const FArrayBox& statein = SBorder[mfi];
        FArrayBox& rhoGradIn = rhoGradFab[mfi];
        computeMaxGrad(statein, rhoGradIn, bx, dx, rhoGradMaxCycle);
        rhoGradMax = fmax(rhoGradMax, rhoGradMaxCycle);
	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
    }
    }
    ParallelDescriptor::ReduceRealMax(rhoGradMax);
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
                Vector<int>  itags;

        const FArrayBox& statein = SBorder[mfi];
        // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
        // So we are going to get a temporary integer array.
        const FArrayBox& rhoGradIn = rhoGradFab[mfi];

        const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    tagfab.get_itags(itags, tilebx);

        Array4<char> tagfabArray = tagfab.array();
            // data pointer and index space

	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

        performTaggingRhoGrad(tagfabArray,  tilebx, statein, rhoGradIn,rhoGradMax, dx);
	    // state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		// 	BL_TO_FORTRAN_3D(S_new[mfi]),
		// 	&tagval, &clearval,
		// 	AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()),
		// 	AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);

	    // Now update the tags in the TagBox.
	    //
	    // tagfab.tags_and_untags(itags, tilebx);
	}
    }
}

void
MultiphaseAMR::computeMaxGrad(const FArrayBox& statein, FArrayBox& rhoGradIn, const Box& bx,
    const Real* dx, Real& rhoGradMax)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(bx);
    Dim3 hi = ubound(bx);
    Real rhoL, rhoC, rhoR, rhoU, rhoD;
    Real dRho_dx, dRho_dy;
    Real temp = 0.0;
    Array4<Real> const& rhoGrad = rhoGradIn.array();

    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {
                rhoL = stateArray(i,j+1,k,0);
                rhoC = stateArray(i+1, j+1, k, 0);
                rhoR = stateArray(i+2, j+1, k, 0);
                rhoU = stateArray(i+1, j+2, k, 0);
                rhoD = stateArray(i+1, j, k, 0);
                dRho_dx = (rhoR-rhoL)/(2*dx[0]);
                dRho_dy = (rhoU-rhoD)/(2*dx[1]);

                rhoGrad(i,j,k,0) = sqrt(dRho_dx*dRho_dx + dRho_dy*dRho_dy);

                if(temp < rhoGrad(i,j,k,0))
                {
                    temp = rhoGrad(i,j,k,0);
                }
            }
        }
    }
    // std::cout << "Break" << " ";

    rhoGradMax = temp;

}

void
MultiphaseAMR::performTaggingRhoGrad(Array4<char> tagfabArray, const Box& tilebx,
    const FArrayBox& statein,  const FArrayBox& rhoGradIn, Real rhoGradMax, const Real* dx)
{
    Array4<Real const> const& stateArray = statein.array();
    Dim3 lo = lbound(tilebx);
    Dim3 hi = ubound(tilebx);
    Real rhoL, rhoC, rhoR, rhoU, rhoD;
    Real dRho_dx, dRho_dy;
    Array4<Real const> const& rhoGrad = rhoGradIn.array();

    for(int i = lo.x; i <= hi.x; i++)
    {
        for(int j = lo.y; j <= hi.y; j++)
        {
            for(int k = lo.z; k <= hi.z; k++)
            {

                if(rhoGrad(i,j,k,0) >= 0.4*rhoGradMax)
                {
                    tagfabArray(i,j,k,0) = TagBox::SET;
                }
                else
                {
                    tagfabArray(i,j,k,0) = TagBox::CLEAR;

                }
            }
        }
    }

}

void
MultiphaseAMR::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }

    // This tutorial code only supports periodic boundaries.
    // if (! gg->isAllPeriodic()) {
	// amrex::Abort("Please set geom.is_periodic = 1 1 1");
    // }

#ifdef AMREX_PARTICLES
    pp.query("do_tracers", do_tracers);
#endif

    //
    // read tagging parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    // get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

void
MultiphaseAMR::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Cons_Type),1.0,0,0,NUM_STATE,geom);

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        amrex::Print() << "MultiphaseAMR::reflux() at level " << level
                       << " : time = " << end << std::endl;
    }
}

void
MultiphaseAMR::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Cons_Type);
}

void
MultiphaseAMR::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    MultiphaseAMR& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);

    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
MultiphaseAMR::init_particles ()
{
  if (do_tracers and level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->do_tiling = true;
      TracerPC->tile_size = IntVect(AMREX_D_DECL(1024000,4,4));

      AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif
