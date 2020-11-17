
subroutine advect(time, lo, hi, &
     &            phiIn , ui_lo, ui_hi, &
     &            phiOut, uo_lo, uo_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flyy, fy_lo, fy_hi, &
     &            dx,dt, nVar) bind(C, name="advect")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_2d

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: nVar
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  double precision, intent(in   ) :: phiIn (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2), 0:nVar-1)
  double precision, intent(inout) :: phiOut(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2), 0:nVar-1)
   double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
  double precision, intent(  out) :: flyy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))

  integer :: i, j, Dimn, n
  integer :: glo(2), ghi(2)
  double precision :: dtdx(2), umax, vmax

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:), pointer, contiguous :: phix_1d, phiy_1d, phix, phiy, slope

  print*, dx(1)
  dtdx = dt/dx

  Dimn = 2
  glo = lo - 1
  ghi = hi + 1

  ! edge states
  call bl_allocate(phix_1d, glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(phiy_1d, glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(phix   , glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(phiy   , glo(1), ghi(1), glo(2), ghi(2))
  ! slope                                                 
  call bl_allocate(slope  , glo(1), ghi(1), glo(2), ghi(2))

  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
!  umax = maxval(abs(vx))
!  vmax = maxval(abs(vy))
!  if ( umax*dt .ge. dx(1) .or. &
!       vmax*dt .ge. dx(2) ) then
!     print *, "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
!     call bl_error("CFL violation. Use smaller adv.cfl.")
!  end if

  ! call a function to compute flux
  call compute_flux_2d(lo, hi, dt, dx, &
                             phiIn,ui_lo,ui_hi, nVar, &
                             flxx, fx_lo, fx_hi, &
                             flyy, fy_lo, fy_hi)

 ! Apply conservative update
  do j = lo(2),hi(2)
	do i = lo(1),hi(1)
		do n = 0, nVar-1
			phiOut(i,j,n) = phiIn(i,j,n) + dtdx(1)*(flxx(i,j) - flxx(i+1, j))
		end do
	end do	
 end do
! Splitting: compute_flux_2d, then apply update, then call flux function again
  ! Final fluxes
!  do    j = lo(2), hi(2)
!     do i = lo(1), hi(1)+1
!        flxx(i,j) = phix(i,j) * vx(i,j)
!     end do
!  end do
!  !
!  do    j = lo(2), hi(2)+1
!     do i = lo(1), hi(1)
!        flxy(i,j) = phiy(i,j) * vy(i,j)
!     end do
!  end do

!  ! Do a conservative update
!  do    j = lo(2),hi(2)
!     do i = lo(1),hi(1)
!        uout(i,j) = uin(i,j) + &
!             ( (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
!             + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) )
!     enddo
!  enddo

!  ! Scale by face area in order to correctly reflx
!  do    j = lo(2), hi(2)
!     do i = lo(1), hi(1)+1
!        flxx(i,j) = flxx(i,j) * ( dt * dx(2))
!     enddo
!  enddo
  
!  ! Scale by face area in order to correctly reflx
!  do    j = lo(2), hi(2)+1 
!     do i = lo(1), hi(1)
!        flxy(i,j) = flxy(i,j) * (dt * dx(1))
!     enddo
!  enddo

!  call bl_deallocate(phix_1d)
!  call bl_deallocate(phiy_1d)
!  call bl_deallocate(phix)
!  call bl_deallocate(phiy)
!  call bl_deallocate(slope)

end subroutine advect

!subroutine computeTimeStep(phi, ph_lo, ph_hi, nVar, &
!	&  dx, gam, dt) bind(C, name="computeTimeStep")
	
!!	use general_functions module, only : Conserv2Prim
	
!	 implicit none
!	 double precision, intent(in) :: dx(2), gam
!     integer, intent(in) :: ph_lo(2), ph_hi(2), nVar
!     double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), 0:nVar-1)
!	 double precision, intent (inout) :: dt
	 
!	 double precision :: Smax_x = 0.0d0
!	 double precision :: Smax_y = 0.0d0
!	 double precision :: phiCell(0:nVar-1)
!	 double precision :: WCell(0:nVar-1)
!	 double precision :: cS
!	 double precision :: CFL = 0.9d0
!	 integer i, j, n

!	 do i = ph_lo(1), ph_hi(1)
!		do j = ph_lo(2), ph_hi(2)
!			do n = 0, nVar-1
!				phiCell = phi(i, j, n)
!			end do
!				call Conserv2Prim(phiCell, gam, WCell, nVar)
!				cS = DSQRT(gam*WCell(3)/WCell(0))
!				Smax_x = DMAX1(Smax_x, cS + DABS(WCell(1)))
!				Smax_y = DMAX1(Smax_y, cS + DABS(WCell(2)))
!		end do
!	end do
!		dt = CFL*DMIN1(dx(1)/Smax_x, dx(2)/Smax_y)		
!	 end subroutine computeTimeStep

!	subroutine Prim2Conserv(W, gam, phi, nVar)
	
!	integer, intent(in) :: nVar
!	double precision, intent(	out) :: phi(0:nVar-1)
!	double precision, intent(in	) :: W(0:nVar-1)
!	double precision, intent(in) :: gam
	
!	phi(0) = W(0)
!	phi(1) = W(0)*W(1)
!	phi(2) = W(0)*W(2)
!	phi(3) = W(3)/(gam - 1.0) + 0.5*W(0)*(W(1)**2.0 + W(2)**2.0)
	
!	end subroutine Prim2Conserv
	
!	subroutine Conserv2Prim(phi, gam, W, nVar)
	
!	integer, intent(in) :: nVar
!	double precision, intent(in	) :: phi(0:nVar-1)
!	double precision, intent(  out) :: W(0:nVar-1)
!	double precision, intent(in) :: gam
	
!	W(0) = phi(0)
!	W(1) = phi(1)/phi(0)
!	W(2) = phi(2)/phi(0)
!	W(3) = (gam - 1.0)*(phi(3) - 0.5*(phi(1)**2.0 + phi(2)**2.0)/phi(0))
	
!	end subroutine Conserv2Prim
