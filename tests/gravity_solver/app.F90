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

    do k=lbound(lgrid%gmgv(iv)%coords_x1,3),ubound(lgrid%gmgv(iv)%coords_x1,3)
     do j=lbound(lgrid%gmgv(iv)%coords_x1,2),ubound(lgrid%gmgv(iv)%coords_x1,2)
      do i=lbound(lgrid%gmgv(iv)%coords_x1,1),ubound(lgrid%gmgv(iv)%coords_x1,1)
       lgrid%gmgv(iv)%coords_x1(i,j,k) = lgrid%gmgv(iv)%nodes(1,i,j,k)
      end do
     end do
    end do

    do k=lbound(lgrid%gmgv(iv)%coords_x2,3),ubound(lgrid%gmgv(iv)%coords_x2,3)
     do j=lbound(lgrid%gmgv(iv)%coords_x2,2),ubound(lgrid%gmgv(iv)%coords_x2,2)
      do i=lbound(lgrid%gmgv(iv)%coords_x2,1),ubound(lgrid%gmgv(iv)%coords_x2,1)
       lgrid%gmgv(iv)%coords_x2(i,j,k) = lgrid%gmgv(iv)%nodes(2,i,j,k)
      end do
     end do
    end do

    do k=lbound(lgrid%gmgv(iv)%coords_x3,3),ubound(lgrid%gmgv(iv)%coords_x3,3)
     do j=lbound(lgrid%gmgv(iv)%coords_x3,2),ubound(lgrid%gmgv(iv)%coords_x3,2)
      do i=lbound(lgrid%gmgv(iv)%coords_x3,1),ubound(lgrid%gmgv(iv)%coords_x3,1)
       lgrid%gmgv(iv)%coords_x3(i,j,k) = lgrid%gmgv(iv)%nodes(3,i,j,k)
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

    do k=lbound(lgrid%gmgv(iv)%vol,3),ubound(lgrid%gmgv(iv)%vol,3)
     do j=lbound(lgrid%gmgv(iv)%vol,2),ubound(lgrid%gmgv(iv)%vol,2)
      do i=lbound(lgrid%gmgv(iv)%vol,1),ubound(lgrid%gmgv(iv)%vol,1)
        lgrid%gmgv(iv)%vol(i,j,k) = &
        (lgrid%gmgv(iv)%coords_x1(i+1,j,k)-lgrid%gmgv(iv)%coords_x1(i,j,k))* &
        (lgrid%gmgv(iv)%coords_x2(i,j+1,k)-lgrid%gmgv(iv)%coords_x2(i,j,k))* &
        (lgrid%gmgv(iv)%coords_x3(i,j,k+1)-lgrid%gmgv(iv)%coords_x3(i,j,k))
      end do
     end do
    end do

   end do

   mgrid%dummy = 0

end subroutine fill_gmg_grids

#endif
