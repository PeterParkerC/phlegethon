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
 real(kind=rp) :: x,y,z,r,r0,rho0,rho
 real(kind=rp) :: x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu
 
 r0 = 0.25_rp
 rho0 = 1.0_rp

 x1l = -0.5_rp
 x1u = +0.5_rp
 x2l = -0.5_rp
 x2u = +0.5_rp
 x3l = -0.5_rp
 x3u = +0.5_rp
 gamma_ad = 5.0_rp/3.0_rp
 mu = 1.0_rp

 call initialize_simulation(mgrid,lgrid,x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu)

 do k=lbound(lgrid%prim,4),ubound(lgrid%prim,4)
  do j=lbound(lgrid%prim,3),ubound(lgrid%prim,3)
   do i=lbound(lgrid%prim,2),ubound(lgrid%prim,2)

     x = lgrid%coords(1,i,j,k)
     y = lgrid%coords(2,i,j,k)
     z = lgrid%coords(3,i,j,k)

     r = sqrt(x**2+y**2+z**2)

     rho = 0.0_rp

     if(r<=r0) then
      rho = rho0*(1.0_rp-r**2/r0**2)**2
     endif

     lgrid%prim(i_rho,i,j,k) = rho
     lgrid%prim(i_vx1,i,j,k) = 0.0_rp
     lgrid%prim(i_vx2,i,j,k) = 0.0_rp
     lgrid%prim(i_vx3,i,j,k) = 0.0_rp
     lgrid%prim(i_p,i,j,k) = 1.0_rp 
 
     lgrid%phi_cc(i,j,k) = 0.0_rp

   end do
  end do
 end do

end subroutine setup_gs

end program test

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
