program test
 use source
 
 type(mpigrid) :: mgrid
 type(locgrid) :: lgrid

 integer :: i,j,k
 real(kind=rp) :: x,y,eta,mach0
 real(kind=rp) :: x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu
 real(kind=rp) :: ye,abar,abar0,dabar

 mach0 = 1e-1_rp
 ye = 0.5_rp
 abar = 0.0_rp
 abar0 = 16.0_rp
 dabar = 0.5_rp

 x1l = 0.0_rp
 x1u = 2.0_rp
 x2l = -0.5_rp
 x2u = 0.5_rp
 x3l = 0.0_rp
 x3u = 1.0_rp
 gamma_ad = 5.0_rp/3.0_rp
 mu = 1.0_rp

 call initialize_simulation(mgrid,lgrid,x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu)

 lgrid%pp_index(1,1)=10
 
 lgrid%pp_index(1,2)=10
 
 lgrid%pp_index(1,3)=1
 

#ifdef USE_MHD

 do k=lbound(lgrid%b_x1,3),ubound(lgrid%b_x1,3)
  do j=lbound(lgrid%b_x1,2),ubound(lgrid%b_x1,2)
   do i=lbound(lgrid%b_x1,1),ubound(lgrid%b_x1,1)
    lgrid%b_x1(i,j,k) = 0.1_rp*mach0
   end do
  end do
 end do

 do k=lbound(lgrid%b_x2,3),ubound(lgrid%b_x2,3)
  do j=lbound(lgrid%b_x2,2),ubound(lgrid%b_x2,2)
   do i=lbound(lgrid%b_x2,1),ubound(lgrid%b_x2,1)
    lgrid%b_x2(i,j,k) = 0.0_rp
   end do
  end do
 end do

#if sdims_make==3

 do k=lbound(lgrid%b_x3,3),ubound(lgrid%b_x3,3)
  do j=lbound(lgrid%b_x3,2),ubound(lgrid%b_x3,2)
   do i=lbound(lgrid%b_x3,1),ubound(lgrid%b_x3,1)
    lgrid%b_x3(i,j,k) = 0.0_rp
   end do
  end do
 end do

#endif

#endif

 do k=lbound(lgrid%prim,4),ubound(lgrid%prim,4)
  do j=lbound(lgrid%prim,3),ubound(lgrid%prim,3)
   do i=lbound(lgrid%prim,2),ubound(lgrid%prim,2)

     x = lgrid%coords(1,i,j,k)
     y = lgrid%coords(2,i,j,k)

     if ((y.gt.(-0.25_rp-1.0_rp/32.0_rp)).and.(y.lt.(-0.25_rp+1.0_rp/32.0_rp))) then
       eta = 0.5_rp*(1.0_rp+sin(16.0_rp*CONST_PI*(y+0.25_rp)))
     else if ((y.ge.(-0.25_rp+1.0_rp/32.0_rp)).and.(y.le.(0.25_rp-1.0_rp/32.0_rp))) then
       eta = 1.0_rp
     else if ((y.gt.(0.25_rp-1.0_rp/32.0_rp)).and.(y.lt.(0.25_rp+1.0_rp/32.0_rp))) then
       eta = 0.5_rp*(1.0_rp+sin(-16.0_rp*CONST_PI*(y-0.25_rp)))
     else
       eta = 0.0
     end if

     lgrid%prim(i_rho,i,j,k) = gamma_ad
     lgrid%prim(i_p,i,j,k) = 1.0_rp

#if nas_make>0
#ifdef ADVECT_YE_IABAR
     abar = abar0*(1.0_rp+dabar*eta)
     lgrid%prim(i_ye,i,j,k) = ye
     lgrid%prim(i_iabar,i,j,k) = 1.0_rp/abar
     mu = abar/(ye*abar+1.0_rp)
#if nas_make>2
     lgrid%prim(i_iabar+1:i_asl,i,j,k) = eta
#endif
#else
     lgrid%prim(i_as1:i_asl,i,j,k) = eta
#endif
#endif

     lgrid%temp(i,j,k) = 1.0_rp/(CONST_RGAS*gamma_ad)*mu

     lgrid%prim(i_vx1,i,j,k) = mach0*(1.0_rp-2.0_rp*eta)
     lgrid%prim(i_vx2,i,j,k) = 0.1_rp*mach0*sin(2.0_rp*CONST_PI*x)
#if sdims_make==3
     lgrid%prim(i_vx3,i,j,k) = 0.0_rp
#endif

   end do
  end do
 end do

 call time_loop(mgrid,lgrid)

 call finalize_simulation(lgrid)
 
end program test

#ifdef USERDEF_OUTPUT

 subroutine extract_userdef_quantities(mgrid,lgrid,iudflush)
 use source        
 type(locgrid), intent(inout) :: lgrid
 type(mpigrid), intent(inout) :: mgrid
 integer, intent(in) :: iudflush


 integer :: i,j,k,ierr
 real(kind=rp) :: Pmax(1),Pmax_comm(1)
 real(kind=rp) :: divb(1),divb_comm(1)
 real(kind=rp) :: ekin(1),ekin_comm(1)
 real(kind=rp) :: dx,abs_b

 !this example extracts: t, max P, location of max P in x and y, mean(nabla B dx / |B|), and Ekin_tot

 Pmax(1) = rp0
 divb(1) = rp0
 ekin(1) = rp0

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)

    Pmax(1) = max(Pmax(1),lgrid%prim(i_p,i,j,k))
    
    abs_b = sqrt(lgrid%b_cc(1,i,j,k)**2+lgrid%b_cc(2,i,j,k)**2)
    dx = rp1/(lgrid%inv_dx1 + lgrid%inv_dx2)
    divb(1) = divb(1) + ( &
    lgrid%inv_dx1*(lgrid%b_x1(i+1,j,k)-lgrid%b_x1(i,j,k)) + &
    lgrid%inv_dx2*(lgrid%b_x2(i,j+1,k)-lgrid%b_x2(i,j,k)) ) / &
    (abs_b/dx)
    
    ekin(1) = ekin(1) + &
    rph*lgrid%prim(i_rho,i,j,k)*( &
    lgrid%prim(i_vx1,i,j,k)**2+lgrid%prim(i_vx2,i,j,k)**2 )

   end do
  end do
 end do

 call mpi_allreduce(Pmax, Pmax_comm, 1, MPI_RP , MPI_MAX, mgrid%comm_cart, ierr)

 call mpi_allreduce(divb, divb_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)

 call mpi_allreduce(ekin, ekin_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)

 divb_comm(1) = divb_comm(1) / real(ddx1*ddx2*ddx3,kind=rp)

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)

    if(abs(lgrid%prim(i_p,i,j,k)/Pmax_comm(1)-rp1)<em11) then
     lgrid%ud_state(iudflush,3) = lgrid%coords(1,i,j,k)
     lgrid%ud_state(iudflush,4) = lgrid%coords(2,i,j,k)
     exit
    endif

   end do
  end do
 end do

 lgrid%ud_state(iudflush,1) = lgrid%time
 lgrid%ud_state(iudflush,2) = Pmax_comm(1)
 lgrid%ud_state(iudflush,5) = divb_comm(1)
 lgrid%ud_state(iudflush,6) = ekin_comm(1)

 mgrid%dummy = 0

 end subroutine extract_userdef_quantities

#endif

