module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d, computeEqnFluxes, computefluxForce

contains

  subroutine computeEqnFluxes(phi, nVar, gam, Dimn, fluxEq)
  
    integer, intent(in) :: nVar, Dimn
    double precision, intent(in) :: gam
    double precision, intent(in   ) :: phi (0:nVar-1)
	double precision, intent(	out) :: fluxEq(0:nVar-1)
  
	double precision :: p
	
	p = (gam-1)*(phi(3) - 0.5*((phi(1)**2 + phi(2)**2))/phi(0))
	if(Dimn.eq. 1) then
		fluxEq(0) = phi(1)
		fluxEq(1) = (phi(1)**2)/phi(0) + p
		fluxEq(2) = phi(1)*(phi(2)/phi(0))
		fluxEq(3) = (phi(1)/phi(0))*(phi(3)+p)
	else if (Dimn.eq. 2) then
		fluxEq(0) = phi(1)
		fluxEq(2) = phi(2)*(phi(1)/phi(0))
		fluxEq(2) = (phi(2)**2)/phi(0) + p
		fluxEq(3) = (phi(1)/phi(0))*(phi(3)+p)
	end if
	
  end subroutine computeEqnFluxes
  
  subroutine computefluxForce(lo, hi, phi, ph_lo, ph_hi, nVar, fluxForce, &
						f_lo, f_hi, dx, dt, gam, Dimn)
  
     double precision, intent(in) :: dt, dx(2)
     double precision, intent(in) :: gam
     integer, intent(in) :: Dimn, nVar
     integer, intent(in) :: ph_lo(2), ph_hi(2)
     integer, intent(in) :: lo(2), hi(2)
	 integer, intent(in) :: f_lo(2), f_hi(2)
	 double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), 0:nVar-1)
	 double precision, intent(	out) :: fluxForce(f_lo(1):f_hi(1), f_lo(2):f_hi(2), 0:nVar-1)

	 integer :: i, j, n
	 double precision :: fL(0:nVar-1)
	 double precision :: fR(0:nVar-1)
	 double precision :: fluxRM(0:nVar-1)
	 double precision :: fluxLF(0:nVar-1)
	 double precision :: phiL(0:nVar-1)
	 double precision :: phiR(0:nVar-1)
	 double precision :: phiMid(0:nVar-1)


      do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1
		do n = 0, nVar - 1
		phiL(n) = phi(i-1, j, n)
		phiR(n) = phi(i, j, n)
		call computeEqnFluxes(phiL, nVar, gam, Dimn, fL)
		call computeEqnFluxes(phiR, nVar, gam, Dimn, fR)
		
		fluxLF(n) = 0.5*(fL(n) - fR(n)) - 0.5*(dx(Dimn)/dt)*(phiR(n) - phiL(n))
		phiMid(n) = 0.5*(phiL(n) - phiR(n)) - 0.5*(dt/dx(Dimn))*(fR(n) - fL(n))
		end do
		
		call computeEqnFluxes(phiMid, nVar, gam, Dimn, fluxRM)
		
		do n = 0, nVar - 1
		fluxForce(i,j,n) = 0.5*(fluxLF(n) + fluxRM(n))
		end do
		
	  end do 
	end do
	
	end subroutine computefluxForce
  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, nVar, &
                             flxx, fx_lo, fx_hi, &
                             flyy, fy_lo, fy_hi)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2) !glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    integer, intent(in) :: nVar
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), 0:nVar-1)
    ! umac, vmac: cell face velocities
    
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2), 0:nVar-1)
    double precision, intent(  out) :: flyy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2), 0:nVar-1)
         
    integer :: i, j
    double precision :: hdtdx(2)
	double precision :: gam
	double precision :: fluxLF(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2), 0:nVar-1)
	double precision :: fluxRM(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2), 0:nVar-1)
	double precision :: fluxForce(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2), 0:nVar-1)

    hdtdx = 0.5*(dt/dx)
	gam = 1.4d0
	
	call computefluxForce(lo, hi, phi, ph_lo, ph_hi, nVar, flxx, &
						fx_lo, fx_hi, dx, dt, gam, 1)
						
	call computefluxForce(lo, hi, phi, ph_lo, ph_hi, nVar, flyy, &
						fy_lo, fy_hi, dx, dt, gam, 2)
						
  end subroutine compute_flux_2d

end module compute_flux_module
