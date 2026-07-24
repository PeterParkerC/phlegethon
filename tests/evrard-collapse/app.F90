program test
 use source
 implicit none

 type(mpigrid) :: mgrid
 type(locgrid) :: lgrid

 call setup_gs(mgrid,lgrid)
 call time_loop(mgrid,lgrid)
 call finalize_simulation(lgrid)

contains
 
subroutine setup_gs(mgrid,lgrid)

 type(mpigrid), intent(inout) :: mgrid
 type(locgrid), intent(inout) :: lgrid

 integer :: i,j,k
 real(kind=rp) :: x,y,z,r,rho
 real(kind=rp) :: x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu
 real(kind=rp) :: uth=0.05_rp,RR=1.0_rp,M=1.0_rp

 x1l = -2.0_rp
 x1u = +2.0_rp
 x2l = -2.0_rp
 x2u = +2.0_rp
 x3l = -2.0_rp
 x3u = +2.0_rp
 gamma_ad = 5.0_rp/3.0_rp
 mu = 1.0_rp

 call initialize_simulation(mgrid,lgrid,x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu)

 lgrid%gs_rho_bg = 0.0_rp

 do k=lbound(lgrid%prim,4),ubound(lgrid%prim,4)
  do j=lbound(lgrid%prim,3),ubound(lgrid%prim,3)
   do i=lbound(lgrid%prim,2),ubound(lgrid%prim,2)

     x = lgrid%coords(1,i,j,k)
     y = lgrid%coords(2,i,j,k)
     z = lgrid%coords(3,i,j,k)

     r = sqrt(x**2+y**2+z**2)

     rho = 1.0e-4_rp

     if(r<=RR) then
      rho = M/(2.0_rp*CONST_PI*RR**2*r)
     endif

     lgrid%prim(i_rho,i,j,k) = rho
     lgrid%prim(i_vx1,i,j,k) = 0.0_rp
     lgrid%prim(i_vx2,i,j,k) = 0.0_rp
     lgrid%prim(i_vx3,i,j,k) = 0.0_rp
     lgrid%prim(i_p,i,j,k) = (gamma_ad-1.0_rp)*rho*uth
     lgrid%phi_cc(i,j,k) = 0.0_rp

   end do
  end do
 end do

end subroutine setup_gs

end program test

subroutine extract_userdef_quantities(mgrid,lgrid,iudflush)
 use source
 type(locgrid), intent(inout) :: lgrid
 type(mpigrid), intent(inout) :: mgrid
 integer, intent(in) :: iudflush

 integer :: i,j,k,ierr
 real(kind=rp) :: etot(1),etot_comm(1)
 real(kind=rp) :: eint(1),eint_comm(1)
 real(kind=rp) :: ekin(1),ekin_comm(1)
 real(kind=rp) :: epot(1),epot_comm(1)
 real(kind=rp) :: a,b,c,d,vol

 eint(1) = rp0
 ekin(1) = rp0
 epot(1) = rp0
 etot(1) = rp0

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)

    vol = 1.0_rp/lgrid%ivol(i,j,k)

    a = lgrid%prim(i_p,i,j,k)/(lgrid%gm-1.0_rp)

    b = 0.5_rp*lgrid%prim(i_rho,i,j,k)*( &
    lgrid%prim(i_vx1,i,j,k)**2 + &
    lgrid%prim(i_vx2,i,j,k)**2 + &
    lgrid%prim(i_vx3,i,j,k)**2 &
    )

    c = 0.5_rp*lgrid%prim(i_rho,i,j,k)*lgrid%phi_cc(i,j,k)
   
    d = a + b + c

    eint(1) = eint(1) + a*vol
    ekin(1) = ekin(1) + b*vol
    epot(1) = epot(1) + c*vol
    etot(1) = etot(1) + d*vol

   end do
  end do
 end do

 call mpi_allreduce(eint, eint_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)
 call mpi_allreduce(ekin, ekin_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)
 call mpi_allreduce(epot, epot_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)
 call mpi_allreduce(etot, etot_comm, 1, MPI_RP , MPI_SUM, mgrid%comm_cart, ierr)

 lgrid%ud_state(iudflush,1) = lgrid%time
 lgrid%ud_state(iudflush,2) = eint_comm(1)
 lgrid%ud_state(iudflush,3) = ekin_comm(1)
 lgrid%ud_state(iudflush,4) = epot_comm(1)
 lgrid%ud_state(iudflush,5) = etot_comm(1)

 mgrid%dummy = 0

end subroutine extract_userdef_quantities

#ifndef GEOMETRY_CARTESIAN_NONUNIFORM

#ifdef GMG_PRECONDITIONER

subroutine fill_gmg_grids(mgrid,lgrid)
   use source
   type(mpigrid), intent(inout) :: mgrid
   type(locgrid), intent(inout) :: lgrid

   integer :: i,j,k,iv,nx1l,nx2l,nx3l
   real(kind=rp) :: xmin,xmax,ymin,ymax,zmin,zmax,dxx,dyy,dzz

   xmin = lgrid%x1l
   xmax = lgrid%x1u
   ymin = lgrid%x2l
   ymax = lgrid%x2u
   zmin = lgrid%x3l
   zmax = lgrid%x3u

   do iv=1,gmg_max_level

    nx1l = int(nx1/2**(iv-1))
    nx2l = int(nx2/2**(iv-1))
    nx3l = int(nx3/2**(iv-1))

    dxx = (xmax-xmin)/real(nx1l,kind=rp)
    dyy = (ymax-ymin)/real(nx2l,kind=rp)
    dzz = (zmax-zmin)/real(nx3l,kind=rp)

    do k = lbound(lgrid%gmgv(iv)%nodes,4), ubound(lgrid%gmgv(iv)%nodes,4)
     do j = lbound(lgrid%gmgv(iv)%nodes,3), ubound(lgrid%gmgv(iv)%nodes,3)
      do i = lbound(lgrid%gmgv(iv)%nodes,2), ubound(lgrid%gmgv(iv)%nodes,2)

       lgrid%gmgv(iv)%nodes(1,i,j,k) = lgrid%x1l + (i-rp1)*dxx
       lgrid%gmgv(iv)%nodes(2,i,j,k) = lgrid%x2l + (j-rp1)*dyy
       lgrid%gmgv(iv)%nodes(3,i,j,k) = lgrid%x3l + (k-rp1)*dzz

      end do
     end do
    end do

    do k=lbound(lgrid%gmgv(iv)%coords,4),ubound(lgrid%gmgv(iv)%coords,4)
     do j=lbound(lgrid%gmgv(iv)%coords,3),ubound(lgrid%gmgv(iv)%coords,3)
      do i=lbound(lgrid%gmgv(iv)%coords,2),ubound(lgrid%gmgv(iv)%coords,2)
       lgrid%gmgv(iv)%coords(1,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(1,i+1,j,k)+lgrid%gmgv(iv)%nodes(1,i,j,k))
       lgrid%gmgv(iv)%coords(2,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(2,i,j+1,k)+lgrid%gmgv(iv)%nodes(2,i,j,k))
       lgrid%gmgv(iv)%coords(3,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(3,i,j,k+1)+lgrid%gmgv(iv)%nodes(3,i,j,k))
      end do
     end do
    end do

   end do

   mgrid%dummy = 0

end subroutine fill_gmg_grids

#endif

#else

subroutine create_geometry(lgrid,mgrid)
   use source
   type(locgrid), intent(inout) :: lgrid
   type(mpigrid), intent(in) :: mgrid

   integer :: i,j,k
   real(kind=rp) :: dxx,dyy,x0,mx,rg,bx,ax,xnorm,xmin,xmax, &
   ay,by,ynorm,ymax,x,y,xx,yy,my,dx,dy,ymin

   real(kind=rp) :: az,bz,znorm,zmax,z,zz,mz,dzz,zmin,dz

   xmin = lgrid%x1l
   xmax = lgrid%x1u
   ymax = lgrid%x2u
   ymin = -ymax
   zmax = lgrid%x3u
   zmin = -zmax

   rg = 2.0_rp

   x0 = 0.0_rp

   mx = 5.0_rp
   my = 5.0_rp

   dxx = 2.0_rp/real(nx1,kind=rp)
   dyy = 2.0_rp/real(nx2,kind=rp)
   dzz = 2.0_rp/real(nx3,kind=rp)
   mz = 5.0_rp

   ax = 1.0_rp/rg
   bx = 1.0_rp - ax
   xnorm = xmax/(ax + bx)

   ay = 1.0_rp/rg
   by = 1.0_rp - ay
   ynorm = ymax/(ay + by)

   az = 1.0_rp/rg
   bz = 1.0_rp - az
   znorm = zmax/(az + bz)

   do k = lbound(lgrid%nodes,4), ubound(lgrid%nodes,4)
    zz = -1.0_rp + (real(k,kind=rp) - 1.0_rp)*dzz
    do j = lbound(lgrid%nodes,3), ubound(lgrid%nodes,3)
     yy = -1.0_rp + (real(j,kind=rp) - 1.0_rp)*dyy
     do i = lbound(lgrid%nodes,2), ubound(lgrid%nodes,2)
      xx = -1.0_rp + (real(i,kind=rp) - 1.0_rp)*dxx
      x = xnorm*(ax*xx + bx*xx**mx)
      y = ynorm*(ay*yy + by*yy**my)
      z = znorm*(az*zz + bz*zz**mz)

      lgrid%nodes(1,i,j,k) = x
      lgrid%nodes(2,i,j,k) = y
      lgrid%nodes(3,i,j,k) = z
     end do
    end do
   end do

   if(mgrid%coords_dd(1)==0) then

     do k=lbound(lgrid%nodes,4),ubound(lgrid%nodes,4)
      do j=lbound(lgrid%nodes,3),ubound(lgrid%nodes,3)
       dx = lgrid%nodes(1,2,j,k)-lgrid%nodes(1,1,j,k)
       do i=1-ngc,0
        lgrid%nodes(1,i,j,k) = xmin + (i-1.0_rp)*dx
       end do
      end do
     end do

   end if

   if(mgrid%coords_dd(1)==mgrid%bricks(1)-1) then

     do k=lbound(lgrid%nodes,4),ubound(lgrid%nodes,4)
      do j=lbound(lgrid%nodes,3),ubound(lgrid%nodes,3)
       dx = lgrid%nodes(1,ubound(lgrid%nodes,2)-ngc,j,k) - &
       lgrid%nodes(1,ubound(lgrid%nodes,2)-ngc-1,j,k)
       do i=ubound(lgrid%nodes,2)-ngc+1,ubound(lgrid%nodes,2)
        lgrid%nodes(1,i,j,k) = xmax + (i-nx1-1.0_rp)*dx
       end do
      end do
     end do

   end if

   if(mgrid%coords_dd(2)==0) then

     do k=lbound(lgrid%nodes,4),ubound(lgrid%nodes,4)
      do i=lbound(lgrid%nodes,2),ubound(lgrid%nodes,2)
       dy = lgrid%nodes(2,i,2,k)-lgrid%nodes(2,i,1,k)
       do j=1-ngc,0
        lgrid%nodes(2,i,j,k) = ymin + (j-1.0_rp)*dy
       end do
      end do
     end do

   end if

   if(mgrid%coords_dd(2)==mgrid%bricks(2)-1) then

     do k=lbound(lgrid%nodes,4),ubound(lgrid%nodes,4)
      do i=lbound(lgrid%nodes,2),ubound(lgrid%nodes,2)
       dy = lgrid%nodes(2,i,ubound(lgrid%nodes,3)-ngc,k) - &
       lgrid%nodes(2,i,ubound(lgrid%nodes,3)-ngc-1,k)
       do j=ubound(lgrid%nodes,3)-ngc+1,ubound(lgrid%nodes,3)
        lgrid%nodes(2,i,j,k) = ymax + (j-nx2-1.0_rp)*dy
       end do
      end do
     end do

   end if

   if(mgrid%coords_dd(3)==0) then

     do j=lbound(lgrid%nodes,3),ubound(lgrid%nodes,3)
      do i=lbound(lgrid%nodes,2),ubound(lgrid%nodes,2)
       dz = lgrid%nodes(3,i,j,2)-lgrid%nodes(3,i,j,1)
       do k=1-ngc,0
        lgrid%nodes(3,i,j,k) = zmin + (k-1.0_rp)*dz
       end do
      end do
     end do

   end if

   if(mgrid%coords_dd(3)==mgrid%bricks(3)-1) then

     do j=lbound(lgrid%nodes,3),ubound(lgrid%nodes,3)
      do i=lbound(lgrid%nodes,2),ubound(lgrid%nodes,2)
       dz = lgrid%nodes(3,i,j,ubound(lgrid%nodes,4)-ngc) - &
       lgrid%nodes(3,i,j,ubound(lgrid%nodes,4)-ngc-1)
       do k=ubound(lgrid%nodes,4)-ngc+1,ubound(lgrid%nodes,4)
        lgrid%nodes(3,i,j,k) = zmax + (k-nx3-1.0_rp)*dz
       end do
      end do
     end do

   end if

   do k=lbound(lgrid%coords_x1,4),ubound(lgrid%coords_x1,4)
    do j=lbound(lgrid%coords_x1,3),ubound(lgrid%coords_x1,3)
     do i=lbound(lgrid%coords_x1,2),ubound(lgrid%coords_x1,2)
      lgrid%coords_x1(1,i,j,k) = lgrid%nodes(1,i,j,k)
      lgrid%coords_x1(2,i,j,k) = 0.5_rp*(lgrid%nodes(2,i,j+1,k)+lgrid%nodes(2,i,j,k))
      lgrid%coords_x1(3,i,j,k) = 0.5_rp*(lgrid%nodes(3,i,j,k+1)+lgrid%nodes(3,i,j,k))
     end do
    end do
   end do

   do k=lbound(lgrid%coords_x2,4),ubound(lgrid%coords_x2,4)
    do j=lbound(lgrid%coords_x2,3),ubound(lgrid%coords_x2,3)
     do i=lbound(lgrid%coords_x2,2),ubound(lgrid%coords_x2,2)
      lgrid%coords_x2(1,i,j,k) = 0.5_rp*(lgrid%nodes(1,i+1,j,k)+lgrid%nodes(1,i,j,k))
      lgrid%coords_x2(2,i,j,k) = lgrid%nodes(2,i,j,k)
      lgrid%coords_x2(3,i,j,k) = 0.5_rp*(lgrid%nodes(3,i,j,k+1)+lgrid%nodes(3,i,j,k))
     end do
    end do
   end do

   do k=lbound(lgrid%coords_x3,4),ubound(lgrid%coords_x3,4)
    do j=lbound(lgrid%coords_x3,3),ubound(lgrid%coords_x3,3)
     do i=lbound(lgrid%coords_x3,2),ubound(lgrid%coords_x3,2)
      lgrid%coords_x3(1,i,j,k) = 0.5_rp*(lgrid%nodes(1,i+1,j,k)+lgrid%nodes(1,i,j,k))
      lgrid%coords_x3(2,i,j,k) = 0.5_rp*(lgrid%nodes(2,i,j+1,k)+lgrid%nodes(2,i,j,k))
      lgrid%coords_x3(3,i,j,k) = lgrid%nodes(3,i,j,k)
     end do
    end do
   end do

   do k=lbound(lgrid%coords,4),ubound(lgrid%coords,4)
    do j=lbound(lgrid%coords,3),ubound(lgrid%coords,3)
     do i=lbound(lgrid%coords,2),ubound(lgrid%coords,2)
      lgrid%coords(1,i,j,k) = 0.5_rp*(lgrid%nodes(1,i+1,j,k)+lgrid%nodes(1,i,j,k))
      lgrid%coords(2,i,j,k) = 0.5_rp*(lgrid%nodes(2,i,j+1,k)+lgrid%nodes(2,i,j,k))
      lgrid%coords(3,i,j,k) = 0.5_rp*(lgrid%nodes(3,i,j,k+1)+lgrid%nodes(3,i,j,k))
     end do
    end do
   end do

end subroutine create_geometry

#ifdef GMG_PRECONDITIONER

subroutine fill_gmg_grids(mgrid,lgrid)
   use source
   type(mpigrid), intent(inout) :: mgrid
   type(locgrid), intent(inout) :: lgrid

   integer :: i,j,k,nx1l,nx2l,nx3l,ngcl
   real(kind=rp) :: dxx,dyy,x0,mx,rg,bx,ax,xnorm,xmin,xmax, &
   ay,by,ynorm,ymax,x,y,xx,yy,my,dx,dy,ymin

   real(kind=rp) :: az,bz,znorm,zmax,z,zz,mz,dzz,zmin,dz

   xmin = lgrid%x1l
   xmax = lgrid%x1u
   ymax = lgrid%x2u
   ymin = lgrid%x2l
   zmax = lgrid%x3u
   zmin = lgrid%x3l

   rg = 2.0_rp

   x0 = 0.0_rp

   mx = 5.0_rp
   my = 5.0_rp
   mz = 5.0_rp

   ax = 1.0_rp/rg
   bx = 1.0_rp - ax
   xnorm = xmax/(ax + bx)

   ay = 1.0_rp/rg
   by = 1.0_rp - ay
   ynorm = ymax/(ay + by)

   az = 1.0_rp/rg
   bz = 1.0_rp - az
   znorm = zmax/(az + bz)

   ngcl = 1

   do iv=1,gmg_max_level

    nx1l = int(nx1/2**(iv-1))
    nx2l = int(nx2/2**(iv-1))
    nx3l = int(nx3/2**(iv-1))

    dxx = 2.0_rp/real(nx1l,kind=rp)
    dyy = 2.0_rp/real(nx2l,kind=rp)
    dzz = 2.0_rp/real(nx3l,kind=rp)

    do k = lbound(lgrid%gmgv(iv)%nodes,4), ubound(lgrid%gmgv(iv)%nodes,4)

     zz = -1.0_rp + (real(k,kind=rp) - 1.0_rp)*dzz

     do j = lbound(lgrid%gmgv(iv)%nodes,3), ubound(lgrid%gmgv(iv)%nodes,3)

      yy = -1.0_rp + (real(j,kind=rp) - 1.0_rp)*dyy

      do i = lbound(lgrid%gmgv(iv)%nodes,2), ubound(lgrid%gmgv(iv)%nodes,2)

       xx = -1.0_rp + (real(i,kind=rp) - 1.0_rp)*dxx
       x = xnorm*(ax*xx + bx*xx**mx)
       y = ynorm*(ay*yy + by*yy**my)
       z = znorm*(az*zz + bz*zz**mz)

       lgrid%gmgv(iv)%nodes(1,i,j,k) = x
       lgrid%gmgv(iv)%nodes(2,i,j,k) = y
       lgrid%gmgv(iv)%nodes(3,i,j,k) = z

      end do
     end do
    end do

    if(mgrid%coords_dd(1)==0) then

     do k=lbound(lgrid%gmgv(iv)%nodes,4),ubound(lgrid%gmgv(iv)%nodes,4)
      do j=lbound(lgrid%gmgv(iv)%nodes,3),ubound(lgrid%gmgv(iv)%nodes,3)
       dx = lgrid%gmgv(iv)%nodes(1,2,j,k)-lgrid%gmgv(iv)%nodes(1,1,j,k)
       do i=1-ngcl,0
        lgrid%gmgv(iv)%nodes(1,i,j,k) = xmin + (i-1.0_rp)*dx
       end do
      end do
      end do

    end if

    if(mgrid%coords_dd(1)==mgrid%bricks(1)-1) then

     do k=lbound(lgrid%gmgv(iv)%nodes,4),ubound(lgrid%gmgv(iv)%nodes,4)
      do j=lbound(lgrid%gmgv(iv)%nodes,3),ubound(lgrid%gmgv(iv)%nodes,3)
       dx = lgrid%gmgv(iv)%nodes(1,ubound(lgrid%gmgv(iv)%nodes,2)-ngcl,j,k) - &
       lgrid%gmgv(iv)%nodes(1,ubound(lgrid%gmgv(iv)%nodes,2)-ngcl-1,j,k)
       do i=ubound(lgrid%gmgv(iv)%nodes,2)-ngcl+1,ubound(lgrid%gmgv(iv)%nodes,2)
        lgrid%gmgv(iv)%nodes(1,i,j,k) = xmax + (i-nx1l-1.0_rp)*dx
       end do
      end do
     end do

    end if

    if(mgrid%coords_dd(2)==0) then

     do k=lbound(lgrid%gmgv(iv)%nodes,4),ubound(lgrid%gmgv(iv)%nodes,4)
      do i=lbound(lgrid%gmgv(iv)%nodes,2),ubound(lgrid%gmgv(iv)%nodes,2)
       dy = lgrid%gmgv(iv)%nodes(2,i,2,k)-lgrid%gmgv(iv)%nodes(2,i,1,k)
       do j=1-ngcl,0
        lgrid%gmgv(iv)%nodes(2,i,j,k) = ymin + (j-1.0_rp)*dy
       end do
      end do
     end do

    end if

    if(mgrid%coords_dd(2)==mgrid%bricks(2)-1) then

     do k=lbound(lgrid%gmgv(iv)%nodes,4),ubound(lgrid%gmgv(iv)%nodes,4)
      do i=lbound(lgrid%gmgv(iv)%nodes,2),ubound(lgrid%gmgv(iv)%nodes,2)
       dy = lgrid%gmgv(iv)%nodes(2,i,ubound(lgrid%gmgv(iv)%nodes,3)-ngcl,k) - &
       lgrid%gmgv(iv)%nodes(2,i,ubound(lgrid%gmgv(iv)%nodes,3)-ngcl-1,k)
       do j=ubound(lgrid%gmgv(iv)%nodes,3)-ngcl+1,ubound(lgrid%gmgv(iv)%nodes,3)
        lgrid%gmgv(iv)%nodes(2,i,j,k) = ymax + (j-nx2l-1.0_rp)*dy
       end do
      end do
     end do

    end if

    if(mgrid%coords_dd(3)==0) then

     do j=lbound(lgrid%gmgv(iv)%nodes,3),ubound(lgrid%gmgv(iv)%nodes,3)
      do i=lbound(lgrid%gmgv(iv)%nodes,2),ubound(lgrid%gmgv(iv)%nodes,2)
       dz = lgrid%gmgv(iv)%nodes(3,i,j,2)-lgrid%gmgv(iv)%nodes(3,i,j,1)
       do k=1-ngcl,0
        lgrid%gmgv(iv)%nodes(3,i,j,k) = zmin + (k-1.0_rp)*dz
       end do
      end do
     end do

    end if

    if(mgrid%coords_dd(3)==mgrid%bricks(3)-1) then

     do j=lbound(lgrid%gmgv(iv)%nodes,3),ubound(lgrid%gmgv(iv)%nodes,3)
      do i=lbound(lgrid%gmgv(iv)%nodes,2),ubound(lgrid%gmgv(iv)%nodes,2)
       dz = lgrid%gmgv(iv)%nodes(3,i,j,ubound(lgrid%gmgv(iv)%nodes,4)-ngcl) - &
       lgrid%gmgv(iv)%nodes(3,i,j,ubound(lgrid%gmgv(iv)%nodes,4)-ngcl-1)
       do k=ubound(lgrid%gmgv(iv)%nodes,4)-ngcl+1,ubound(lgrid%gmgv(iv)%nodes,4)
        lgrid%gmgv(iv)%nodes(3,i,j,k) = zmax + (k-nx3l-1.0_rp)*dz
       end do
      end do
     end do

    end if

    do k=lbound(lgrid%gmgv(iv)%coords,4),ubound(lgrid%gmgv(iv)%coords,4)
     do j=lbound(lgrid%gmgv(iv)%coords,3),ubound(lgrid%gmgv(iv)%coords,3)
      do i=lbound(lgrid%gmgv(iv)%coords,2),ubound(lgrid%gmgv(iv)%coords,2)
       lgrid%gmgv(iv)%coords(1,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(1,i+1,j,k)+lgrid%gmgv(iv)%nodes(1,i,j,k))
       lgrid%gmgv(iv)%coords(2,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(2,i,j+1,k)+lgrid%gmgv(iv)%nodes(2,i,j,k))
       lgrid%gmgv(iv)%coords(3,i,j,k) = 0.5_rp*(lgrid%gmgv(iv)%nodes(3,i,j,k+1)+lgrid%gmgv(iv)%nodes(3,i,j,k))
      end do
     end do
    end do

   end do

   mgrid%dummy = 0

end subroutine fill_gmg_grids

#endif

#endif
