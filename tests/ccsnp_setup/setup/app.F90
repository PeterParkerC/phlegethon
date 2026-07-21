program ccsn_progenitor
 use source
 implicit none

 type model_t
  real(kind=rp), allocatable :: r(:),grav(:),gpot(:),p(:),rho(:),T(:),Xs(:,:)
 end type model_t

 type(mpigrid) :: mgrid
 type(locgrid) :: lgrid

 call setup(mgrid,lgrid)

 call time_loop(mgrid,lgrid)

 call finalize_simulation(lgrid)

contains

subroutine setup(mgrid,lgrid)
 type(mpigrid),intent(inout) :: mgrid
 type(locgrid),intent(inout) :: lgrid

 type(model_t) :: model
 integer :: i,j,k,iv
 real(kind=rp) :: x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad=5.0_rp/3.0_rp,mu=1.0_rp
 real(kind=rp) :: rho,T,ye,abar,zbar,eint,cv,c,p,inv_abar,sound,grav,gpot

 real(kind=rp), dimension(1:nspecies) :: Xs
 real(kind=rp) :: x,y,z=0.0_rp,r,theta,phi,rmin,rmax,rnd

#ifdef USE_MHD
 real(kind=rp) :: Rcy,acy,x0,y0,Aphi,A0
 real(kind=rp), allocatable, dimension(:,:,:) :: Ax,Ay
#endif

 Xs(:) = 0.0_rp

 call load_model('ccsn_progenitor.in',model)

 rmin = 0.6e9_rp 
 rmax = 0.6e10_rp 

#ifdef USE_MHD
 Rcy = 0.55_rp*rmax
 acy = 0.1_rp*rmax
 A0 = 1e4_rp*acy
#endif

 x1l = -rmax
 x1u = +rmax
 x2l = -rmax
 x2u = +rmax
 x3l = -rmax
 x3u = +rmax

 call initialize_simulation(mgrid,lgrid,x1l,x1u,x2l,x2u,x3l,x3u,gamma_ad,mu)

#ifdef SAVE_PLANES
 lgrid%planes_x3_index(1) = int(nx3/2)
#endif

#ifdef USE_MHD

 allocate(Ax(mgrid%i1(1):mgrid%i2(1)+1,mgrid%i1(2):mgrid%i2(2)+1,mgrid%i1(3):mgrid%i2(3)+1))
 allocate(Ay(mgrid%i1(1):mgrid%i2(1)+1,mgrid%i1(2):mgrid%i2(2)+1,mgrid%i1(3):mgrid%i2(3)+1))

 do k=mgrid%i1(3),mgrid%i2(3)+1
  do j=mgrid%i1(2),mgrid%i2(2)+1
   do i=mgrid%i1(1),mgrid%i2(1)+1

     x = lgrid%nodes(1,i,j,k)
     y = lgrid%nodes(2,i,j,k)
     z = lgrid%nodes(3,i,j,k)
     phi = atan2(y,x)

     x0 = Rcy*cos(phi)
     y0 = Rcy*sin(phi)

     r = sqrt((x-x0)**2+(y-y0)**2+z**2)/acy

     Ax(i,j,k) = rp0
     Ay(i,j,k) = rp0

     if(r<=rp1) then

      Aphi = A0*cos(0.5_rp*CONST_PI*r)**2
      Ax(i,j,k) = -Aphi*sin(phi)
      Ay(i,j,k) =  Aphi*cos(phi)

     end if

   end do
  end do
 end do

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)+1

     lgrid%b_x1(i,j,k) = -(Ay(i,j,k+1)-Ay(i,j,k))*lgrid%inv_dx1

   end do
  end do
 end do

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)+1
   do i=mgrid%i1(1),mgrid%i2(1)

     lgrid%b_x2(i,j,k) = (Ax(i,j,k+1)-Ax(i,j,k))*lgrid%inv_dx1

   end do
  end do
 end do 
               
 do k=mgrid%i1(3),mgrid%i2(3)+1
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)

     lgrid%b_x3(i,j,k) = (Ay(i+1,j,k)-Ay(i,j,k)-Ax(i,j+1,k)+Ax(i,j,k))*lgrid%inv_dx1

   end do
  end do
 end do

 deallocate(Ax)
 deallocate(Ay)

#endif

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)+1

     x = lgrid%coords_x1(1,i,j,k)
     y = lgrid%coords_x1(2,i,j,k)
     z = lgrid%coords_x1(3,i,j,k)

     r = sqrt(x**2+y**2+z**2)

     call interpolate_model(model,r,grav,gpot,p,rho,T,Xs)

     ye = 0.0_rp
     inv_abar = 0.0_rp
     do iv=1,nspecies
      inv_abar = inv_abar + Xs(iv)/lgrid%A(iv)
      ye = ye + Xs(iv)*lgrid%Z(iv)/lgrid%A(iv)
     end do
     abar = 1.0_rp/inv_abar
     zbar = ye*abar

     call helm_rhoT_given_full(rho,T,abar,zbar,p,eint,sound,cv)

#ifdef EVOLVE_ETOT
     lgrid%phi_x1(i,j,k) = gpot
#endif

#ifdef USE_WB
     lgrid%eq_prim_x1(ieq_rho,i,j,k) = rho
     lgrid%eq_prim_x1(ieq_p,i,j,k) = p
     lgrid%eq_prim_x1(ieq_T,i,j,k) = T
#ifdef USE_FASTEOS
     eint = rho*eint
     lgrid%eq_gammae_x1(i,j,k) = p/eint+1.0_rp
#endif
#endif

   end do
  end do
 end do

 do k=mgrid%i1(3),mgrid%i2(3)
  do j=mgrid%i1(2),mgrid%i2(2)+1
   do i=mgrid%i1(1),mgrid%i2(1)

     x = lgrid%coords_x2(1,i,j,k)
     y = lgrid%coords_x2(2,i,j,k)
     z = lgrid%coords_x2(3,i,j,k)

     r = sqrt(x**2+y**2+z**2)

     call interpolate_model(model,r,grav,gpot,p,rho,T,Xs)

     ye = 0.0_rp
     inv_abar = 0.0_rp
     do iv=1,nspecies
      inv_abar = inv_abar + Xs(iv)/lgrid%A(iv)
      ye = ye + Xs(iv)*lgrid%Z(iv)/lgrid%A(iv)
     end do
     abar = 1.0_rp/inv_abar
     zbar = ye*abar

     call helm_rhoT_given_full(rho,T,abar,zbar,p,eint,sound,cv)

#ifdef EVOLVE_ETOT
     lgrid%phi_x2(i,j,k) = gpot
#endif

#ifdef USE_WB
     lgrid%eq_prim_x2(ieq_rho,i,j,k) = rho
     lgrid%eq_prim_x2(ieq_p,i,j,k) = p
     lgrid%eq_prim_x2(ieq_T,i,j,k) = T
#ifdef USE_FASTEOS
     eint = rho*eint
     lgrid%eq_gammae_x2(i,j,k) = p/eint+1.0_rp
#endif
#endif

   end do
  end do
 end do

 do k=mgrid%i1(3),mgrid%i2(3)+1
  do j=mgrid%i1(2),mgrid%i2(2)
   do i=mgrid%i1(1),mgrid%i2(1)

     x = lgrid%coords_x3(1,i,j,k)
     y = lgrid%coords_x3(2,i,j,k)
     z = lgrid%coords_x3(3,i,j,k)

     r = sqrt(x**2+y**2+z**2)

     call interpolate_model(model,r,grav,gpot,p,rho,T,Xs)

     ye = 0.0_rp
     inv_abar = 0.0_rp
     do iv=1,nspecies
      inv_abar = inv_abar + Xs(iv)/lgrid%A(iv)
      ye = ye + Xs(iv)*lgrid%Z(iv)/lgrid%A(iv)
     end do
     abar = 1.0_rp/inv_abar
     zbar = ye*abar

     call helm_rhoT_given_full(rho,T,abar,zbar,p,eint,sound,cv)

#ifdef EVOLVE_ETOT
     lgrid%phi_x3(i,j,k) = gpot
#endif

#ifdef USE_WB
     lgrid%eq_prim_x3(ieq_rho,i,j,k) = rho
     lgrid%eq_prim_x3(ieq_p,i,j,k) = p
     lgrid%eq_prim_x3(ieq_T,i,j,k) = T
#ifdef USE_FASTEOS
     eint = rho*eint
     lgrid%eq_gammae_x3(i,j,k) = p/eint+1.0_rp
#endif
#endif

   end do
  end do
 end do

 do k=lbound(lgrid%prim,4),ubound(lgrid%prim,4)
  do j=lbound(lgrid%prim,3),ubound(lgrid%prim,3)
   do i=lbound(lgrid%prim,2),ubound(lgrid%prim,2)

     x = lgrid%coords(1,i,j,k)
     y = lgrid%coords(2,i,j,k)
     z = lgrid%coords(3,i,j,k)
     r = sqrt(x**2+y**2+z**2)

     theta = acos(z/r)
     phi = atan2(y,x)

     call interpolate_model(model,r,grav,gpot,p,rho,T,Xs)

     call helm_rhoT_given_full(rho,T,abar,zbar,p,eint,c,cv)

     lgrid%prim(i_rho,i,j,k) = rho

     call random_number(rnd)
     lgrid%prim(i_p,i,j,k) = p*(1.0_rp+1.0e-4_rp*(2.0_rp*rnd-1.0_rp))

     lgrid%temp(i,j,k) = T

     lgrid%prim(i_vx1,i,j,k) = 0.0_rp
     lgrid%prim(i_vx2,i,j,k) = 0.0_rp
     lgrid%prim(i_vx3,i,j,k) = 0.0_rp

     ye = 0.0_rp
     inv_abar = 0.0_rp
     do iv=1,nspecies
      inv_abar = inv_abar + Xs(iv)/lgrid%A(iv)
      ye = ye + Xs(iv)*lgrid%Z(iv)/lgrid%A(iv)
     end do 
     abar = 1.0_rp/inv_abar
     zbar = ye*abar

     lgrid%prim(i_as1:,i,j,k) = 0.0_rp

     do iv=1,nspecies
      lgrid%prim(i_as1+iv-1,i,j,k) = Xs(iv)
     end do

     lgrid%grav(1,i,j,k) = grav*cos(phi)*sin(theta)
     lgrid%grav(2,i,j,k) = grav*sin(phi)*sin(theta)
     lgrid%grav(3,i,j,k) = grav*cos(theta)

#ifdef EVOLVE_ETOT
     lgrid%phi_cc(i,j,k) = gpot
#endif

#ifdef USE_WB
     lgrid%eq_prim_cc(ieq_rho,i,j,k) = rho
     lgrid%eq_prim_cc(ieq_p,i,j,k) = p
     lgrid%eq_prim_cc(ieq_T,i,j,k) = T
#ifdef USE_FASTEOS
     eint = rho*eint
     lgrid%eq_gammae_cc(i,j,k) = p/eint+1.0_rp
#endif
#endif

#ifdef USE_INTERNAL_BOUNDARIES
    if ((r>=rmin) .and. (r<=rmax)) then
     lgrid%is_solid(i,j,k) = 0
    else
     lgrid%is_solid(i,j,k) = 1
    endif
#endif

   end do
  end do
 end do

 call deallocate_model(model)

end subroutine setup

subroutine load_model(filename, m)
  character(len=*), intent(in) :: filename
  type(model_t), intent(inout) :: m
  integer, parameter :: fp=23
  integer :: i, len

  open(unit=fp, file=filename)
  read(fp,*) len
  allocate(m%r(len),m%grav(len),m%gpot(len),m%p(len),m%rho(len), &
   m%T(len),m%Xs(len,1:nspecies))

  do i=1, len
    read(fp, *) m%r(i),m%grav(i),m%gpot(i),m%p(i),m%rho(i),m%T(i),m%Xs(i,:)
  end do

  close(fp)

end subroutine load_model

subroutine interpolate_model(m, r, grav, gpot, p, rho, T, Xs)
  type(model_t), intent(in) :: m
  real(kind=rp), intent(in) :: r
  real(kind=rp), intent(out) :: grav, gpot, p, rho, T
  real(kind=rp), dimension(1:nspecies), intent(out) :: Xs

  real(kind=rp) :: fac, rr
  integer :: i, iv

  rr = r

  if (rr < m%r(lbound(m%r,1))) then
    rr = m%r(lbound(m%r,1))
  end if

  if (rr > m%r(ubound(m%r,1))) then
    rr = m%r(ubound(m%r,1))
  end if

  i = minloc(abs(m%r(:) - rr),1)
  if (m%r(i) > rr) i = i - 1
  i = min(max(lbound(m%r,1),i),ubound(m%r,1)-1)

  fac = (rr - m%r(i)) / (m%r(i+1) - m%r(i))
  grav = fac * (m%grav(i+1) - m%grav(i)) + m%grav(i)
  gpot = fac * (m%gpot(i+1) - m%gpot(i)) + m%gpot(i)
  p = fac * (m%p(i+1) - m%p(i)) + m%p(i)
  rho = fac * (m%rho(i+1) - m%rho(i)) + m%rho(i)
  T = fac * (m%T(i+1) - m%T(i)) + m%T(i)
  do iv=1,nspecies
   Xs(iv) = fac * (m%Xs(i+1,iv) - m%Xs(i,iv)) + m%Xs(i,iv)
  end do

end subroutine interpolate_model

subroutine deallocate_model(m)
  type(model_t), intent(inout) :: m
  deallocate(m%r, m%grav, m%gpot, m%p, m%rho, m%T, m%Xs)
end subroutine deallocate_model

end program ccsn_progenitor

#ifdef USE_NUCLEAR_NETWORK 
subroutine extract_network_information(lgrid) 
  use source 
  type(locgrid), intent(inout) :: lgrid 

  lgrid%A(1)=1.0_rp 
  lgrid%A(2)=4.0_rp 
  lgrid%A(3)=12.0_rp 
  lgrid%A(4)=16.0_rp 
  lgrid%A(5)=20.0_rp 
  lgrid%A(6)=23.0_rp 
  lgrid%A(7)=24.0_rp 
  lgrid%A(8)=28.0_rp 
  lgrid%A(9)=31.0_rp 
  lgrid%A(10)=32.0_rp 
  lgrid%A(11)=36.0_rp 
  lgrid%A(12)=40.0_rp 

  lgrid%Z(1)=1.0_rp 
  lgrid%Z(2)=2.0_rp 
  lgrid%Z(3)=6.0_rp 
  lgrid%Z(4)=8.0_rp 
  lgrid%Z(5)=10.0_rp 
  lgrid%Z(6)=11.0_rp 
  lgrid%Z(7)=12.0_rp 
  lgrid%Z(8)=14.0_rp 
  lgrid%Z(9)=15.0_rp 
  lgrid%Z(10)=16.0_rp 
  lgrid%Z(11)=18.0_rp 
  lgrid%Z(12)=20.0_rp 

  lgrid%name_species(1)='p' 
  lgrid%name_species(2)='he4' 
  lgrid%name_species(3)='c12' 
  lgrid%name_species(4)='o16' 
  lgrid%name_species(5)='ne20' 
  lgrid%name_species(6)='na23' 
  lgrid%name_species(7)='mg24' 
  lgrid%name_species(8)='si28' 
  lgrid%name_species(9)='p31' 
  lgrid%name_species(10)='s32' 
  lgrid%name_species(11)='ar36' 
  lgrid%name_species(12)='ca40' 

  lgrid%name_reacs(1)='o16-->he4+c12' 
  lgrid%name_reacs(2)='ne20-->he4+o16' 
  lgrid%name_reacs(3)='mg24-->p+na23' 
  lgrid%name_reacs(4)='mg24-->he4+ne20' 
  lgrid%name_reacs(5)='si28-->he4+mg24' 
  lgrid%name_reacs(6)='s32-->p+p31' 
  lgrid%name_reacs(7)='s32-->he4+si28' 
  lgrid%name_reacs(8)='ar36-->he4+s32' 
  lgrid%name_reacs(9)='ca40-->he4+ar36' 
  lgrid%name_reacs(10)='c12-->he4+he4+he4' 
  lgrid%name_reacs(11)='he4+c12-->o16' 
  lgrid%name_reacs(12)='he4+o16-->ne20' 
  lgrid%name_reacs(13)='he4+ne20-->mg24' 
  lgrid%name_reacs(14)='p+na23-->mg24' 
  lgrid%name_reacs(15)='he4+mg24-->si28' 
  lgrid%name_reacs(16)='he4+si28-->s32' 
  lgrid%name_reacs(17)='p+p31-->s32' 
  lgrid%name_reacs(18)='he4+s32-->ar36' 
  lgrid%name_reacs(19)='he4+ar36-->ca40' 
  lgrid%name_reacs(20)='c12+c12-->p+na23' 
  lgrid%name_reacs(21)='c12+c12-->he4+ne20' 
  lgrid%name_reacs(22)='c12+o16-->he4+mg24' 
  lgrid%name_reacs(23)='o16+o16-->p+p31' 
  lgrid%name_reacs(24)='o16+o16-->he4+si28' 
  lgrid%name_reacs(25)='he4+ne20-->p+na23' 
  lgrid%name_reacs(26)='he4+ne20-->c12+c12' 
  lgrid%name_reacs(27)='c12+ne20-->p+p31' 
  lgrid%name_reacs(28)='c12+ne20-->he4+si28' 
  lgrid%name_reacs(29)='p+na23-->he4+ne20' 
  lgrid%name_reacs(30)='p+na23-->c12+c12' 
  lgrid%name_reacs(31)='he4+mg24-->c12+o16' 
  lgrid%name_reacs(32)='he4+si28-->p+p31' 
  lgrid%name_reacs(33)='he4+si28-->c12+ne20' 
  lgrid%name_reacs(34)='he4+si28-->o16+o16' 
  lgrid%name_reacs(35)='p+p31-->he4+si28' 
  lgrid%name_reacs(36)='p+p31-->c12+ne20' 
  lgrid%name_reacs(37)='p+p31-->o16+o16' 
  lgrid%name_reacs(38)='he4+he4+he4-->c12' 

#ifdef USE_LMP_WEAK_RATES 
 
  allocate(lgrid%weak_table(1:0,1:13,1:11)) 
  allocate(lgrid%weak_neu(1:0,1:13,1:11)) 
  allocate(lgrid%neu_rates(1:0)) 

#endif 

#ifdef PARTITION_FUNCTIONS_FOR_REVERSE_RATES 
 
  lgrid%temp_part(1)=0.1000_rp 
  lgrid%temp_part(2)=0.1500_rp 
  lgrid%temp_part(3)=0.2000_rp 
  lgrid%temp_part(4)=0.3000_rp 
  lgrid%temp_part(5)=0.4000_rp 
  lgrid%temp_part(6)=0.5000_rp 
  lgrid%temp_part(7)=0.6000_rp 
  lgrid%temp_part(8)=0.7000_rp 
  lgrid%temp_part(9)=0.8000_rp 
  lgrid%temp_part(10)=0.9000_rp 
  lgrid%temp_part(11)=1.0000_rp 
  lgrid%temp_part(12)=1.5000_rp 
  lgrid%temp_part(13)=2.0000_rp 
  lgrid%temp_part(14)=2.5000_rp 
  lgrid%temp_part(15)=3.0000_rp 
  lgrid%temp_part(16)=3.5000_rp 
  lgrid%temp_part(17)=4.0000_rp 
  lgrid%temp_part(18)=4.5000_rp 
  lgrid%temp_part(19)=5.0000_rp 
  lgrid%temp_part(20)=6.0000_rp 
  lgrid%temp_part(21)=7.0000_rp 
  lgrid%temp_part(22)=8.0000_rp 
  lgrid%temp_part(23)=9.0000_rp 
  lgrid%temp_part(24)=10.0000_rp 

  allocate(lgrid%part(1:nspecies,1:24)) 

  lgrid%part(1,1)=0.000e+00_rp 
  lgrid%part(1,2)=0.000e+00_rp 
  lgrid%part(1,3)=0.000e+00_rp 
  lgrid%part(1,4)=0.000e+00_rp 
  lgrid%part(1,5)=0.000e+00_rp 
  lgrid%part(1,6)=0.000e+00_rp 
  lgrid%part(1,7)=0.000e+00_rp 
  lgrid%part(1,8)=0.000e+00_rp 
  lgrid%part(1,9)=0.000e+00_rp 
  lgrid%part(1,10)=0.000e+00_rp 
  lgrid%part(1,11)=0.000e+00_rp 
  lgrid%part(1,12)=0.000e+00_rp 
  lgrid%part(1,13)=0.000e+00_rp 
  lgrid%part(1,14)=0.000e+00_rp 
  lgrid%part(1,15)=0.000e+00_rp 
  lgrid%part(1,16)=0.000e+00_rp 
  lgrid%part(1,17)=0.000e+00_rp 
  lgrid%part(1,18)=0.000e+00_rp 
  lgrid%part(1,19)=0.000e+00_rp 
  lgrid%part(1,20)=0.000e+00_rp 
  lgrid%part(1,21)=0.000e+00_rp 
  lgrid%part(1,22)=0.000e+00_rp 
  lgrid%part(1,23)=0.000e+00_rp 
  lgrid%part(1,24)=0.000e+00_rp 

  lgrid%part(2,1)=0.000e+00_rp 
  lgrid%part(2,2)=0.000e+00_rp 
  lgrid%part(2,3)=0.000e+00_rp 
  lgrid%part(2,4)=0.000e+00_rp 
  lgrid%part(2,5)=0.000e+00_rp 
  lgrid%part(2,6)=0.000e+00_rp 
  lgrid%part(2,7)=0.000e+00_rp 
  lgrid%part(2,8)=0.000e+00_rp 
  lgrid%part(2,9)=0.000e+00_rp 
  lgrid%part(2,10)=0.000e+00_rp 
  lgrid%part(2,11)=0.000e+00_rp 
  lgrid%part(2,12)=0.000e+00_rp 
  lgrid%part(2,13)=0.000e+00_rp 
  lgrid%part(2,14)=0.000e+00_rp 
  lgrid%part(2,15)=0.000e+00_rp 
  lgrid%part(2,16)=0.000e+00_rp 
  lgrid%part(2,17)=0.000e+00_rp 
  lgrid%part(2,18)=0.000e+00_rp 
  lgrid%part(2,19)=0.000e+00_rp 
  lgrid%part(2,20)=0.000e+00_rp 
  lgrid%part(2,21)=0.000e+00_rp 
  lgrid%part(2,22)=0.000e+00_rp 
  lgrid%part(2,23)=0.000e+00_rp 
  lgrid%part(2,24)=0.000e+00_rp 

  lgrid%part(3,1)=0.000e+00_rp 
  lgrid%part(3,2)=0.000e+00_rp 
  lgrid%part(3,3)=0.000e+00_rp 
  lgrid%part(3,4)=0.000e+00_rp 
  lgrid%part(3,5)=0.000e+00_rp 
  lgrid%part(3,6)=0.000e+00_rp 
  lgrid%part(3,7)=0.000e+00_rp 
  lgrid%part(3,8)=0.000e+00_rp 
  lgrid%part(3,9)=0.000e+00_rp 
  lgrid%part(3,10)=0.000e+00_rp 
  lgrid%part(3,11)=0.000e+00_rp 
  lgrid%part(3,12)=0.000e+00_rp 
  lgrid%part(3,13)=0.000e+00_rp 
  lgrid%part(3,14)=0.000e+00_rp 
  lgrid%part(3,15)=0.000e+00_rp 
  lgrid%part(3,16)=0.000e+00_rp 
  lgrid%part(3,17)=0.000e+00_rp 
  lgrid%part(3,18)=0.000e+00_rp 
  lgrid%part(3,19)=0.000e+00_rp 
  lgrid%part(3,20)=0.000e+00_rp 
  lgrid%part(3,21)=0.000e+00_rp 
  lgrid%part(3,22)=0.000e+00_rp 
  lgrid%part(3,23)=0.000e+00_rp 
  lgrid%part(3,24)=0.000e+00_rp 

  lgrid%part(4,1)=0.000e+00_rp 
  lgrid%part(4,2)=0.000e+00_rp 
  lgrid%part(4,3)=0.000e+00_rp 
  lgrid%part(4,4)=0.000e+00_rp 
  lgrid%part(4,5)=0.000e+00_rp 
  lgrid%part(4,6)=0.000e+00_rp 
  lgrid%part(4,7)=0.000e+00_rp 
  lgrid%part(4,8)=0.000e+00_rp 
  lgrid%part(4,9)=0.000e+00_rp 
  lgrid%part(4,10)=0.000e+00_rp 
  lgrid%part(4,11)=0.000e+00_rp 
  lgrid%part(4,12)=0.000e+00_rp 
  lgrid%part(4,13)=0.000e+00_rp 
  lgrid%part(4,14)=0.000e+00_rp 
  lgrid%part(4,15)=0.000e+00_rp 
  lgrid%part(4,16)=0.000e+00_rp 
  lgrid%part(4,17)=0.000e+00_rp 
  lgrid%part(4,18)=0.000e+00_rp 
  lgrid%part(4,19)=0.000e+00_rp 
  lgrid%part(4,20)=0.000e+00_rp 
  lgrid%part(4,21)=0.000e+00_rp 
  lgrid%part(4,22)=0.000e+00_rp 
  lgrid%part(4,23)=0.000e+00_rp 
  lgrid%part(4,24)=0.000e+00_rp 

  lgrid%part(5,1)=0.000e+00_rp 
  lgrid%part(5,2)=0.000e+00_rp 
  lgrid%part(5,3)=0.000e+00_rp 
  lgrid%part(5,4)=0.000e+00_rp 
  lgrid%part(5,5)=0.000e+00_rp 
  lgrid%part(5,6)=0.000e+00_rp 
  lgrid%part(5,7)=0.000e+00_rp 
  lgrid%part(5,8)=0.000e+00_rp 
  lgrid%part(5,9)=0.000e+00_rp 
  lgrid%part(5,10)=0.000e+00_rp 
  lgrid%part(5,11)=0.000e+00_rp 
  lgrid%part(5,12)=1.600e-05_rp 
  lgrid%part(5,13)=3.819e-04_rp 
  lgrid%part(5,14)=2.541e-03_rp 
  lgrid%part(5,15)=8.963e-03_rp 
  lgrid%part(5,16)=2.197e-02_rp 
  lgrid%part(5,17)=4.282e-02_rp 
  lgrid%part(5,18)=7.155e-02_rp 
  lgrid%part(5,19)=1.073e-01_rp 
  lgrid%part(5,20)=1.949e-01_rp 
  lgrid%part(5,21)=2.952e-01_rp 
  lgrid%part(5,22)=4.014e-01_rp 
  lgrid%part(5,23)=5.094e-01_rp 
  lgrid%part(5,24)=6.179e-01_rp 

  lgrid%part(6,1)=0.000e+00_rp 
  lgrid%part(6,2)=0.000e+00_rp 
  lgrid%part(6,3)=0.000e+00_rp 
  lgrid%part(6,4)=0.000e+00_rp 
  lgrid%part(6,5)=4.000e-06_rp 
  lgrid%part(6,6)=5.500e-05_rp 
  lgrid%part(6,7)=3.020e-04_rp 
  lgrid%part(6,8)=1.018e-03_rp 
  lgrid%part(6,9)=2.533e-03_rp 
  lgrid%part(6,10)=5.140e-03_rp 
  lgrid%part(6,11)=9.048e-03_rp 
  lgrid%part(6,12)=4.865e-02_rp 
  lgrid%part(6,13)=1.104e-01_rp 
  lgrid%part(6,14)=1.779e-01_rp 
  lgrid%part(6,15)=2.424e-01_rp 
  lgrid%part(6,16)=3.012e-01_rp 
  lgrid%part(6,17)=3.544e-01_rp 
  lgrid%part(6,18)=4.030e-01_rp 
  lgrid%part(6,19)=4.484e-01_rp 
  lgrid%part(6,20)=5.330e-01_rp 
  lgrid%part(6,21)=6.141e-01_rp 
  lgrid%part(6,22)=6.947e-01_rp 
  lgrid%part(6,23)=7.763e-01_rp 
  lgrid%part(6,24)=8.597e-01_rp 

  lgrid%part(7,1)=0.000e+00_rp 
  lgrid%part(7,2)=0.000e+00_rp 
  lgrid%part(7,3)=0.000e+00_rp 
  lgrid%part(7,4)=0.000e+00_rp 
  lgrid%part(7,5)=0.000e+00_rp 
  lgrid%part(7,6)=0.000e+00_rp 
  lgrid%part(7,7)=0.000e+00_rp 
  lgrid%part(7,8)=0.000e+00_rp 
  lgrid%part(7,9)=0.000e+00_rp 
  lgrid%part(7,10)=0.000e+00_rp 
  lgrid%part(7,11)=1.000e-06_rp 
  lgrid%part(7,12)=1.260e-04_rp 
  lgrid%part(7,13)=1.776e-03_rp 
  lgrid%part(7,14)=8.665e-03_rp 
  lgrid%part(7,15)=2.479e-02_rp 
  lgrid%part(7,16)=5.210e-02_rp 
  lgrid%part(7,17)=9.017e-02_rp 
  lgrid%part(7,18)=1.370e-01_rp 
  lgrid%part(7,19)=1.902e-01_rp 
  lgrid%part(7,20)=3.068e-01_rp 
  lgrid%part(7,21)=4.272e-01_rp 
  lgrid%part(7,22)=5.457e-01_rp 
  lgrid%part(7,23)=6.609e-01_rp 
  lgrid%part(7,24)=7.730e-01_rp 

  lgrid%part(8,1)=0.000e+00_rp 
  lgrid%part(8,2)=0.000e+00_rp 
  lgrid%part(8,3)=0.000e+00_rp 
  lgrid%part(8,4)=0.000e+00_rp 
  lgrid%part(8,5)=0.000e+00_rp 
  lgrid%part(8,6)=0.000e+00_rp 
  lgrid%part(8,7)=0.000e+00_rp 
  lgrid%part(8,8)=0.000e+00_rp 
  lgrid%part(8,9)=0.000e+00_rp 
  lgrid%part(8,10)=0.000e+00_rp 
  lgrid%part(8,11)=0.000e+00_rp 
  lgrid%part(8,12)=5.000e-06_rp 
  lgrid%part(8,13)=1.640e-04_rp 
  lgrid%part(8,14)=1.295e-03_rp 
  lgrid%part(8,15)=5.119e-03_rp 
  lgrid%part(8,16)=1.362e-02_rp 
  lgrid%part(8,17)=2.828e-02_rp 
  lgrid%part(8,18)=4.968e-02_rp 
  lgrid%part(8,19)=7.761e-02_rp 
  lgrid%part(8,20)=1.497e-01_rp 
  lgrid%part(8,21)=2.366e-01_rp 
  lgrid%part(8,22)=3.313e-01_rp 
  lgrid%part(8,23)=4.297e-01_rp 
  lgrid%part(8,24)=5.302e-01_rp 

  lgrid%part(9,1)=0.000e+00_rp 
  lgrid%part(9,2)=0.000e+00_rp 
  lgrid%part(9,3)=0.000e+00_rp 
  lgrid%part(9,4)=0.000e+00_rp 
  lgrid%part(9,5)=0.000e+00_rp 
  lgrid%part(9,6)=0.000e+00_rp 
  lgrid%part(9,7)=0.000e+00_rp 
  lgrid%part(9,8)=0.000e+00_rp 
  lgrid%part(9,9)=0.000e+00_rp 
  lgrid%part(9,10)=0.000e+00_rp 
  lgrid%part(9,11)=1.000e-06_rp 
  lgrid%part(9,12)=1.110e-04_rp 
  lgrid%part(9,13)=1.295e-03_rp 
  lgrid%part(9,14)=5.683e-03_rp 
  lgrid%part(9,15)=1.536e-02_rp 
  lgrid%part(9,16)=3.152e-02_rp 
  lgrid%part(9,17)=5.451e-02_rp 
  lgrid%part(9,18)=8.416e-02_rp 
  lgrid%part(9,19)=1.201e-01_rp 
  lgrid%part(9,20)=2.092e-01_rp 
  lgrid%part(9,21)=3.198e-01_rp 
  lgrid%part(9,22)=4.515e-01_rp 
  lgrid%part(9,23)=6.053e-01_rp 
  lgrid%part(9,24)=7.825e-01_rp 

  lgrid%part(10,1)=0.000e+00_rp 
  lgrid%part(10,2)=0.000e+00_rp 
  lgrid%part(10,3)=0.000e+00_rp 
  lgrid%part(10,4)=0.000e+00_rp 
  lgrid%part(10,5)=0.000e+00_rp 
  lgrid%part(10,6)=0.000e+00_rp 
  lgrid%part(10,7)=0.000e+00_rp 
  lgrid%part(10,8)=0.000e+00_rp 
  lgrid%part(10,9)=0.000e+00_rp 
  lgrid%part(10,10)=0.000e+00_rp 
  lgrid%part(10,11)=0.000e+00_rp 
  lgrid%part(10,12)=0.000e+00_rp 
  lgrid%part(10,13)=1.200e-05_rp 
  lgrid%part(10,14)=1.600e-04_rp 
  lgrid%part(10,15)=8.966e-04_rp 
  lgrid%part(10,16)=3.078e-03_rp 
  lgrid%part(10,17)=7.779e-03_rp 
  lgrid%part(10,18)=1.603e-02_rp 
  lgrid%part(10,19)=2.868e-02_rp 
  lgrid%part(10,20)=6.934e-02_rp 
  lgrid%part(10,21)=1.324e-01_rp 
  lgrid%part(10,22)=2.188e-01_rp 
  lgrid%part(10,23)=3.292e-01_rp 
  lgrid%part(10,24)=4.643e-01_rp 

  lgrid%part(11,1)=0.000e+00_rp 
  lgrid%part(11,2)=0.000e+00_rp 
  lgrid%part(11,3)=0.000e+00_rp 
  lgrid%part(11,4)=0.000e+00_rp 
  lgrid%part(11,5)=0.000e+00_rp 
  lgrid%part(11,6)=0.000e+00_rp 
  lgrid%part(11,7)=0.000e+00_rp 
  lgrid%part(11,8)=0.000e+00_rp 
  lgrid%part(11,9)=0.000e+00_rp 
  lgrid%part(11,10)=0.000e+00_rp 
  lgrid%part(11,11)=0.000e+00_rp 
  lgrid%part(11,12)=1.000e-06_rp 
  lgrid%part(11,13)=5.400e-05_rp 
  lgrid%part(11,14)=5.329e-04_rp 
  lgrid%part(11,15)=2.446e-03_rp 
  lgrid%part(11,16)=7.263e-03_rp 
  lgrid%part(11,17)=1.643e-02_rp 
  lgrid%part(11,18)=3.103e-02_rp 
  lgrid%part(11,19)=5.170e-02_rp 
  lgrid%part(11,20)=1.123e-01_rp 
  lgrid%part(11,21)=1.990e-01_rp 
  lgrid%part(11,22)=3.123e-01_rp 
  lgrid%part(11,23)=4.527e-01_rp 
  lgrid%part(11,24)=6.210e-01_rp 

  lgrid%part(12,1)=0.000e+00_rp 
  lgrid%part(12,2)=0.000e+00_rp 
  lgrid%part(12,3)=0.000e+00_rp 
  lgrid%part(12,4)=0.000e+00_rp 
  lgrid%part(12,5)=0.000e+00_rp 
  lgrid%part(12,6)=0.000e+00_rp 
  lgrid%part(12,7)=0.000e+00_rp 
  lgrid%part(12,8)=0.000e+00_rp 
  lgrid%part(12,9)=0.000e+00_rp 
  lgrid%part(12,10)=0.000e+00_rp 
  lgrid%part(12,11)=0.000e+00_rp 
  lgrid%part(12,12)=0.000e+00_rp 
  lgrid%part(12,13)=0.000e+00_rp 
  lgrid%part(12,14)=0.000e+00_rp 
  lgrid%part(12,15)=8.000e-06_rp 
  lgrid%part(12,16)=6.000e-05_rp 
  lgrid%part(12,17)=2.860e-04_rp 
  lgrid%part(12,18)=9.795e-04_rp 
  lgrid%part(12,19)=2.655e-03_rp 
  lgrid%part(12,20)=1.226e-02_rp 
  lgrid%part(12,21)=3.788e-02_rp 
  lgrid%part(12,22)=9.068e-02_rp 
  lgrid%part(12,23)=1.821e-01_rp 
  lgrid%part(12,24)=3.217e-01_rp 

#endif 
 
end subroutine extract_network_information 

subroutine compute_jina_rates(T9,R) 
  use source 
  real(kind=rp), intent(in) :: T9 
  real(kind=rp), dimension(1:nreacs), intent(inout) :: R 

  real(kind=rp) :: logT9,cf,e1,e2,e3,e4,e5,tmp,tmpp 

  logT9 = log(T9) 
  cf = rp1/T9**fvthirds 
  e1 = T9**tthirds 
  e2 = e1*e1 
  e3 = e2*e1
  e4 = e3*e1 
  e5 = e4*e1 
  e1 = e1*cf 
  e2 = e2*cf 
  e3 = e3*cf
  e4 = e4*cf 
  e5 = e5*cf 
  tmp = rp0 
  tmpp = rp0 

  tmp=rp0 
  tmpp =9.431310e+01_rp & 
  -8.450300e+01_rp*e1 & 
  +5.891280e+01_rp*e2 & 
  -1.482730e+02_rp*e3 & 
  +9.083240e+00_rp*e4 & 
  -5.410410e-01_rp*e5 & 
  +7.185540e+01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.792950e+02_rp & 
  -8.495150e+01_rp*e1 & 
  +1.034110e+02_rp*e2 & 
  -4.205670e+02_rp*e3 & 
  +6.408740e+01_rp*e4 & 
  -1.246240e+01_rp*e5 & 
  +1.388030e+02_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(1)=tmp 
 
  tmp=rp0 
  tmpp =4.866040e+01_rp & 
  -5.488750e+01_rp*e1 & 
  -3.972620e+01_rp*e2 & 
  -2.107990e-01_rp*e3 & 
  +4.428790e-01_rp*e4 & 
  -7.977530e-02_rp*e5 & 
  +8.333330e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.864310e+01_rp & 
  -6.524600e+01_rp*e1
  tmp=tmp+exp(tmpp) 
  tmpp =3.426580e+01_rp & 
  -6.765180e+01_rp*e1 & 
  -3.659250e+00_rp*e3 & 
  +7.142240e-01_rp*e4 & 
  -1.075080e-03_rp*e5
  tmp=tmp+exp(tmpp) 
  R(2)=tmp 
 
  tmp=rp0 
  tmpp =1.684969e+01_rp & 
  -1.372372e+02_rp*e1
  tmp=tmp+exp(tmpp) 
  tmpp =3.458400e+01_rp & 
  -1.389617e+02_rp*e1 & 
  -2.328100e-01_rp*e2 & 
  -6.444300e-01_rp*e3 & 
  +1.398780e+00_rp*e4 & 
  -1.817700e-01_rp*e5
  tmp=tmp+exp(tmpp) 
  tmpp =4.250278e+01_rp & 
  -1.356938e+02_rp*e1 & 
  -2.089222e+01_rp*e2 & 
  +3.273860e+00_rp*e3 & 
  +2.382420e+00_rp*e4 & 
  -9.557700e-01_rp*e5 & 
  +4.265300e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(3)=tmp 
 
  tmp=rp0 
  tmpp =4.932440e+01_rp & 
  -1.081140e+02_rp*e1 & 
  -4.625250e+01_rp*e2 & 
  +5.589010e+00_rp*e3 & 
  +7.618430e+00_rp*e4 & 
  -3.683000e+00_rp*e5 & 
  +8.333330e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.680170e+01_rp & 
  -1.173340e+02_rp*e1
  tmp=tmp+exp(tmpp) 
  tmpp =1.602030e+01_rp & 
  -1.208950e+02_rp*e1 & 
  +1.692290e+01_rp*e3 & 
  -2.573250e+00_rp*e4 & 
  +2.089970e-01_rp*e5
  tmp=tmp+exp(tmpp) 
  R(4)=tmp 
 
  tmp=rp0 
  tmpp =-5.863647e+03_rp & 
  -1.227700e+03_rp*e1 & 
  +1.603558e+04_rp*e2 & 
  -9.292824e+03_rp*e3 & 
  +2.539336e+02_rp*e4 & 
  -8.091800e+00_rp*e5 & 
  +7.132553e+03_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =3.556075e+01_rp & 
  -2.028227e+02_rp*e1 & 
  +1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(5)=tmp 
 
  tmp=rp0 
  tmpp =4.361090e+01_rp & 
  -1.028600e+02_rp*e1 & 
  -2.532780e+01_rp*e2 & 
  +6.493100e+00_rp*e3 & 
  -9.275130e+00_rp*e4 & 
  -6.104390e-01_rp*e5 & 
  +8.333330e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.168290e+01_rp & 
  -1.051190e+02_rp*e1
  tmp=tmp+exp(tmpp) 
  tmpp =2.517290e+01_rp & 
  -1.066370e+02_rp*e1 & 
  +8.093410e+00_rp*e3 & 
  -6.159710e-01_rp*e4 & 
  +3.115900e-02_rp*e5
  tmp=tmp+exp(tmpp) 
  R(6)=tmp 
 
  tmp=7.281300e+01_rp & 
  -8.062600e+01_rp*e1 & 
  -5.948960e+01_rp*e2 & 
  +4.472050e+00_rp*e3 & 
  -4.789890e+00_rp*e4 & 
  +5.572010e-01_rp*e5 & 
  +8.333330e-01_rp*logT9
  R(7)=exp(tmp) 

  tmp=7.381640e+01_rp & 
  -7.706270e+01_rp*e1 & 
  -6.537090e+01_rp*e2 & 
  +5.682940e+00_rp*e3 & 
  -5.003880e+00_rp*e4 & 
  +5.714070e-01_rp*e5 & 
  +8.333330e-01_rp*logT9
  R(8)=exp(tmp) 

  tmp=7.728260e+01_rp & 
  -8.169160e+01_rp*e1 & 
  -7.100460e+01_rp*e2 & 
  +4.065600e+00_rp*e3 & 
  -5.265090e+00_rp*e4 & 
  +6.835460e-01_rp*e5 & 
  +8.333330e-01_rp*logT9
  R(9)=exp(tmp) 

  tmp=2.239400e+01_rp & 
  -8.854930e+01_rp*e1 & 
  -1.349000e+01_rp*e2 & 
  +2.142590e+01_rp*e3 & 
  -1.347690e+00_rp*e4 & 
  +8.798160e-02_rp*e5 & 
  -1.016530e+01_rp*logT9
  R(10)=exp(tmp) 

  tmp=rp0 
  tmpp =6.965260e+01_rp & 
  -1.392540e+00_rp*e1 & 
  +5.891280e+01_rp*e2 & 
  -1.482730e+02_rp*e3 & 
  +9.083240e+00_rp*e4 & 
  -5.410410e-01_rp*e5 & 
  +7.035540e+01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.546340e+02_rp & 
  -1.840970e+00_rp*e1 & 
  +1.034110e+02_rp*e2 & 
  -4.205670e+02_rp*e3 & 
  +6.408740e+01_rp*e4 & 
  -1.246240e+01_rp*e5 & 
  +1.373030e+02_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(11)=tmp 
 
  tmp=rp0 
  tmpp =2.390300e+01_rp & 
  -3.972620e+01_rp*e2 & 
  -2.107990e-01_rp*e3 & 
  +4.428790e-01_rp*e4 & 
  -7.977530e-02_rp*e5 & 
  -6.666670e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =3.885710e+00_rp & 
  -1.035850e+01_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =9.508480e+00_rp & 
  -1.276430e+01_rp*e1 & 
  -3.659250e+00_rp*e3 & 
  +7.142240e-01_rp*e4 & 
  -1.075080e-03_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(12)=tmp 
 
  tmp=rp0 
  tmpp =2.450580e+01_rp & 
  -4.625250e+01_rp*e2 & 
  +5.589010e+00_rp*e3 & 
  +7.618430e+00_rp*e4 & 
  -3.683000e+00_rp*e5 & 
  -6.666670e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-3.870550e+01_rp & 
  -2.506050e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =1.983070e+00_rp & 
  -9.220260e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-8.798270e+00_rp & 
  -1.278090e+01_rp*e1 & 
  +1.692290e+01_rp*e3 & 
  -2.573250e+00_rp*e4 & 
  +2.089970e-01_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(13)=tmp 
 
  tmp=rp0 
  tmpp =-8.178447e+00_rp & 
  -1.543412e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-6.013458e+01_rp & 
  -5.838050e+00_rp*e1 & 
  +4.388380e+01_rp*e2 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =9.555867e+00_rp & 
  -3.267900e+00_rp*e1 & 
  -2.328100e-01_rp*e2 & 
  -6.444300e-01_rp*e3 & 
  +1.398780e+00_rp*e4 & 
  -1.817700e-01_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =1.747465e+01_rp & 
  -2.089222e+01_rp*e2 & 
  +3.273860e+00_rp*e3 & 
  +2.382420e+00_rp*e4 & 
  -9.557700e-01_rp*e5 & 
  -1.073470e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(14)=tmp 
 
  tmp=rp0 
  tmpp =-5.888508e+03_rp & 
  -1.111834e+03_rp*e1 & 
  +1.603558e+04_rp*e2 & 
  -9.292824e+03_rp*e3 & 
  +2.539336e+02_rp*e4 & 
  -8.091800e+00_rp*e5 & 
  +7.131053e+03_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =1.070000e+01_rp & 
  -8.695672e+01_rp*e1
  tmp=tmp+exp(tmpp) 
  R(15)=tmp 
 
  tmp=4.792120e+01_rp & 
  -5.948960e+01_rp*e2 & 
  +4.472050e+00_rp*e3 & 
  -4.789890e+00_rp*e4 & 
  +5.572010e-01_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(16)=exp(tmp) 

  tmp=rp0 
  tmpp =1.925960e+01_rp & 
  -2.532780e+01_rp*e2 & 
  +6.493100e+00_rp*e3 & 
  -9.275130e+00_rp*e4 & 
  -6.104390e-01_rp*e5 & 
  -6.666670e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-2.668390e+00_rp & 
  -2.259580e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =8.215560e-01_rp & 
  -3.777040e+00_rp*e1 & 
  +8.093410e+00_rp*e3 & 
  -6.159710e-01_rp*e4 & 
  +3.115900e-02_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(17)=tmp 
 
  tmp=4.890100e+01_rp & 
  -6.537090e+01_rp*e2 & 
  +5.682940e+00_rp*e3 & 
  -5.003880e+00_rp*e4 & 
  +5.714070e-01_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(18)=exp(tmp) 

  tmp=5.234860e+01_rp & 
  -7.100460e+01_rp*e2 & 
  +4.065600e+00_rp*e3 & 
  -5.265090e+00_rp*e4 & 
  +6.835460e-01_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(19)=exp(tmp) 

  tmp=6.096490e+01_rp & 
  -8.416500e+01_rp*e2 & 
  -1.419100e+00_rp*e3 & 
  -1.146190e-01_rp*e4 & 
  -7.030700e-02_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(20)=exp(tmp) 

  tmp=6.128630e+01_rp & 
  -8.416500e+01_rp*e2 & 
  -1.566270e+00_rp*e3 & 
  -7.360840e-02_rp*e4 & 
  -7.279700e-02_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(21)=exp(tmp) 

  tmp=4.853410e+01_rp & 
  +3.720400e-01_rp*e1 & 
  -1.334130e+02_rp*e2 & 
  +5.015720e+01_rp*e3 & 
  -3.159870e+00_rp*e4 & 
  +1.782510e-02_rp*e5 & 
  -2.370270e+01_rp*logT9
  R(22)=exp(tmp) 

  tmp=8.526280e+01_rp & 
  +2.234530e-01_rp*e1 & 
  -1.458440e+02_rp*e2 & 
  +8.726120e+00_rp*e3 & 
  -5.540350e-01_rp*e4 & 
  -1.375620e-01_rp*e5 & 
  -6.888070e+00_rp*logT9
  R(23)=exp(tmp) 

  tmp=9.724350e+01_rp & 
  -2.685140e-01_rp*e1 & 
  -1.193240e+02_rp*e2 & 
  -3.224970e+01_rp*e3 & 
  +1.462140e+00_rp*e4 & 
  -2.008930e-01_rp*e5 & 
  +1.321480e+01_rp*logT9
  R(24)=exp(tmp) 

  tmp=rp0 
  tmpp =1.918520e+01_rp & 
  -2.757380e+01_rp*e1 & 
  -2.000240e+01_rp*e2 & 
  +1.159880e+01_rp*e3 & 
  -1.373980e+00_rp*e4 & 
  -1.000000e+00_rp*e5 & 
  -6.666670e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =2.274720e-01_rp & 
  -2.943480e+01_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-6.377720e+00_rp & 
  -2.988960e+01_rp*e1 & 
  +1.972970e+01_rp*e3 & 
  -2.209870e+00_rp*e4 & 
  +1.533740e-01_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(25)=tmp 
 
  tmp=6.147480e+01_rp & 
  -5.362670e+01_rp*e1 & 
  -8.416500e+01_rp*e2 & 
  -1.566270e+00_rp*e3 & 
  -7.360840e-02_rp*e4 & 
  -7.279700e-02_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(26)=exp(tmp) 

  tmp=-2.681360e+02_rp & 
  -3.876240e+01_rp*e1 & 
  +3.611540e+02_rp*e2 & 
  -9.264300e+01_rp*e3 & 
  -9.987380e+00_rp*e4 & 
  +8.927370e-01_rp*e5 & 
  +1.610420e+02_rp*logT9
  R(27)=exp(tmp) 

  tmp=-3.089050e+02_rp & 
  -4.721750e+01_rp*e1 & 
  +5.141970e+02_rp*e2 & 
  -2.008960e+02_rp*e3 & 
  -6.427130e+00_rp*e4 & 
  +7.582560e-01_rp*e5 & 
  +2.363590e+02_rp*logT9
  R(28)=exp(tmp) 

  tmp=rp0 
  tmpp =1.897560e+01_rp & 
  -2.000240e+01_rp*e2 & 
  +1.159880e+01_rp*e3 & 
  -1.373980e+00_rp*e4 & 
  -1.000000e+00_rp*e5 & 
  -6.666670e-01_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =1.782950e-02_rp & 
  -1.861030e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-6.587360e+00_rp & 
  -2.315770e+00_rp*e1 & 
  +1.972970e+01_rp*e3 & 
  -2.209870e+00_rp*e4 & 
  +1.533740e-01_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(29)=tmp 
 
  tmp=6.094380e+01_rp & 
  -2.601840e+01_rp*e1 & 
  -8.416500e+01_rp*e2 & 
  -1.419100e+00_rp*e3 & 
  -1.146190e-01_rp*e4 & 
  -7.030700e-02_rp*e5 & 
  -6.666670e-01_rp*logT9
  R(30)=exp(tmp) 

  tmp=4.957380e+01_rp & 
  -7.820200e+01_rp*e1 & 
  -1.334130e+02_rp*e2 & 
  +5.015720e+01_rp*e3 & 
  -3.159870e+00_rp*e4 & 
  +1.782510e-02_rp*e5 & 
  -2.370270e+01_rp*logT9
  R(31)=exp(tmp) 

  tmp=rp0 
  tmpp =-1.143350e+01_rp & 
  -2.566060e+01_rp*e1 & 
  +2.152100e+01_rp*e3 & 
  -1.903550e+00_rp*e4 & 
  +9.272400e-02_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-1.345950e+01_rp & 
  -2.411200e+01_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(32)=tmp 
 
  tmp=-3.077620e+02_rp & 
  -1.867220e+02_rp*e1 & 
  +5.141970e+02_rp*e2 & 
  -2.008960e+02_rp*e3 & 
  -6.427130e+00_rp*e4 & 
  +7.582560e-01_rp*e5 & 
  +2.363590e+02_rp*logT9
  R(33)=exp(tmp) 

  tmp=9.779040e+01_rp & 
  -1.115950e+02_rp*e1 & 
  -1.193240e+02_rp*e2 & 
  -3.224970e+01_rp*e3 & 
  +1.462140e+00_rp*e4 & 
  -2.008930e-01_rp*e5 & 
  +1.321480e+01_rp*logT9
  R(34)=exp(tmp) 

  tmp=rp0 
  tmpp =-1.291900e+01_rp & 
  -1.877160e+00_rp*e1 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  tmpp =-1.089300e+01_rp & 
  -3.425750e+00_rp*e1 & 
  +2.152100e+01_rp*e3 & 
  -1.903550e+00_rp*e4 & 
  +9.272400e-02_rp*e5 & 
  -1.500000e+00_rp*logT9
  tmp=tmp+exp(tmpp) 
  R(35)=tmp 
 
  tmp=-2.664520e+02_rp & 
  -1.560190e+02_rp*e1 & 
  +3.611540e+02_rp*e2 & 
  -9.264300e+01_rp*e3 & 
  -9.987380e+00_rp*e4 & 
  +8.927370e-01_rp*e5 & 
  +1.610420e+02_rp*logT9
  R(36)=exp(tmp) 

  tmp=8.635010e+01_rp & 
  -8.887970e+01_rp*e1 & 
  -1.458440e+02_rp*e2 & 
  +8.726120e+00_rp*e3 & 
  -5.540350e-01_rp*e4 & 
  -1.375620e-01_rp*e5 & 
  -6.888070e+00_rp*logT9
  R(37)=exp(tmp) 

  tmp=-2.435050e+01_rp & 
  -4.126560e+00_rp*e1 & 
  -1.349000e+01_rp*e2 & 
  +2.142590e+01_rp*e3 & 
  -1.347690e+00_rp*e4 & 
  +8.798160e-02_rp*e5 & 
  -1.316530e+01_rp*logT9
  R(38)=exp(tmp) 

end subroutine compute_jina_rates 

#ifdef USE_LMP_WEAK_RATES 
subroutine compute_weak_rates(rhoye,T9,dt,weak_table,weak_neu,neu_rates,R) 
  use source 
  real(kind=rp), intent(in) :: rhoye, T9, dt 
  real(kind=rp), dimension(1:0,1:13,1:11), intent(in) :: weak_table 
  real(kind=rp), dimension(1:0,1:13,1:11), intent(in) :: weak_neu 
  real(kind=rp), dimension(1:0), intent(inout) :: neu_rates 
  real(kind=rp), dimension(1:nreacs), intent(inout) :: R 


  real(kind=rp) :: logrhoye,tmp,u,v,omu,omv,omu_omv,u_omv,v_omu,u_v 
  real(kind=rp) :: lam_neu,E_neu_avg 
  integer :: i,idx_logrhoye,idx_T9
 
  logrhoye = log10(rhoye) 
  u = rp0 
  v = rp0 
  idx_T9 = 0 
  idx_logrhoye = 0 

  do i=1,12 
   if((T9>=weak_T9(i)) .and. (T9<weak_T9(i+1))) then 
    idx_T9 = i 
    tmp = weak_T9(i) 
    u = (T9-tmp)/(weak_T9(i+1)-tmp) 
    exit 
   end if 
  end do 

  do i=1,10 
   if((logrhoye>=weak_logrhoye(i)) .and. (logrhoye<weak_logrhoye(i+1))) then 
    idx_logrhoye = i 
    tmp = weak_logrhoye(i) 
    v = (logrhoye-tmp)/(weak_logrhoye(i+1)-tmp) 
    exit 
   end if 
  end do 

  omu=rp1-u 
  omv=rp1-v 
  omu_omv=omu*omv 
  u_omv=u*omv 
  v_omu=v*omu 
  u_v=u*v 

end subroutine compute_weak_rates 

subroutine compute_weak_neuloss(rho,Y,neu_rates,dedt) 
  use source 
  real(kind=rp), intent(in) :: rho 
  real(kind=rp), dimension(1:nspecies), intent(in) :: Y 
  real(kind=rp), dimension(1:0), intent(in) :: neu_rates 
  real(kind=rp), dimension(1:nreacs), intent(inout) :: dedt 

  real(kind=rp) :: fac 

  fac = rho*CONST_NAV_MEV_TO_ERG 

end subroutine compute_weak_neuloss 
#endif 

#ifdef PARTITION_FUNCTIONS_FOR_REVERSE_RATES 
subroutine use_partition_functions(T9,temp_part,part,R) 
  use source 
  real(kind=rp), intent(in) :: T9 
  real(kind=rp), intent(in) :: temp_part(1:24) 
  real(kind=rp), intent(in) :: part(1:nspecies,1:24) 
  real(kind=rp), dimension(1:nreacs), intent(inout) :: R 

  real(kind=rp) :: fac,tmp,tmpp 
  integer :: idx_temp,i 
  real(kind=rp) :: & 
  part1, & 
  part2, & 
  part3, & 
  part4, & 
  part5, & 
  part6, & 
  part7, & 
  part8, & 
  part9, & 
  part10, & 
  part11, & 
  part12 
  idx_temp = 1 
  fac = rp1 
  tmp = rp1 
  tmpp = rp1 

  do i=1,23 
   if((T9>=temp_part(i)) .and. (T9<temp_part(i+1))) then 
    idx_temp = i 
    tmp = temp_part(i) 
    fac = (T9-tmp)/(temp_part(i+1)-tmp) 
    exit 
   end if 
  end do 

  part1=rp1 
  part2=rp1 
  part3=rp1 
  part4=rp1 
  tmp=part(5,idx_temp) 
  tmpp=tmp+fac*(part(5,idx_temp+1)-tmp) 
  part5=exp(tmpp) 
  tmp=part(6,idx_temp) 
  tmpp=tmp+fac*(part(6,idx_temp+1)-tmp) 
  part6=exp(tmpp) 
  tmp=part(7,idx_temp) 
  tmpp=tmp+fac*(part(7,idx_temp+1)-tmp) 
  part7=exp(tmpp) 
  tmp=part(8,idx_temp) 
  tmpp=tmp+fac*(part(8,idx_temp+1)-tmp) 
  part8=exp(tmpp) 
  tmp=part(9,idx_temp) 
  tmpp=tmp+fac*(part(9,idx_temp+1)-tmp) 
  part9=exp(tmpp) 
  tmp=part(10,idx_temp) 
  tmpp=tmp+fac*(part(10,idx_temp+1)-tmp) 
  part10=exp(tmpp) 
  tmp=part(11,idx_temp) 
  tmpp=tmp+fac*(part(11,idx_temp+1)-tmp) 
  part11=exp(tmpp) 
  tmp=part(12,idx_temp) 
  tmpp=tmp+fac*(part(12,idx_temp+1)-tmp) 
  part12=exp(tmpp) 

  tmp=part4
  R(1)=R(1)*part2*part3/tmp

  tmp=part5
  R(2)=R(2)*part2*part4/tmp

  tmp=part7
  R(3)=R(3)*part1*part6/tmp

  tmp=part7
  R(4)=R(4)*part2*part5/tmp

  tmp=part8
  R(5)=R(5)*part2*part7/tmp

  tmp=part10
  R(6)=R(6)*part1*part9/tmp

  tmp=part10
  R(7)=R(7)*part2*part8/tmp

  tmp=part11
  R(8)=R(8)*part2*part10/tmp

  tmp=part12
  R(9)=R(9)*part2*part11/tmp

  tmp=part3
  R(10)=R(10)*part2*part2*part2/tmp

  tmp=part2*part5
  R(25)=R(25)*part1*part6/tmp

  tmp=part2*part5
  R(26)=R(26)*part3*part3/tmp

  tmp=part1*part6
  R(30)=R(30)*part3*part3/tmp

  tmp=part2*part7
  R(31)=R(31)*part3*part4/tmp

  tmp=part2*part8
  R(32)=R(32)*part1*part9/tmp

  tmp=part2*part8
  R(33)=R(33)*part3*part5/tmp

  tmp=part2*part8
  R(34)=R(34)*part4*part4/tmp

  tmp=part1*part9
  R(36)=R(36)*part3*part5/tmp

  tmp=part1*part9
  R(37)=R(37)*part4*part4/tmp

end subroutine use_partition_functions 
#endif 

#ifdef USE_ELECTRON_SCREENING 
subroutine screen_rates(rho,T,Y,rates) 
  use source 
  real(kind=rp), intent(in) :: T 
  real(kind=rp), intent(in) :: rho 
  real(kind=rp), intent(in) :: Y(1:nspecies) 
  real(kind=rp), intent(inout) :: rates(1:nreacs) 

  real(kind=rp) :: ye,iabar,abar,zbar,Zk,Yk,fw 
  real(kind=rp) :: rhoye,gammap,gammapo4,lngammap,iT,iT3 
  real(kind=rp) :: gammaeff,Z1,Z2,A1,A2,Ahat,Zprod,Zsum 
  real(kind=rp) :: C,zt,jt,tau,b,b3,b4,b5,b6,Hs,Hw,H 
  real(kind=rp) :: Zsum_5_12,Z1_5_12,Z2_5_12 
  real(kind=rp) :: Zsum_5_3,Z1_5_3,Z2_5_3 

  Hw = rp0 
  A1 = rp0 
  A2 = rp0 
  A1 = rp0 
  Ahat = rp0 
  Z1_5_12 = rp0 
  Z2_5_12 = rp0 
  Zsum_5_12 = rp0 
  Z1_5_3 = rp0 
  Z2_5_3 = rp0 
  Zsum_5_3 = rp0 
  zt = rp0 
  jt = rp0 
  C = rp0 
  tau = rp0 
  b = rp0 
  b3 = rp0 
  b4 = rp0 
  b5 = rp0 
  b6 = rp0 
  Hs = rp0 

  ye = rp0 
  iabar = rp0 
  fw = rp0 
  Yk = Y(1) 
  Zk = 1.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(2) 
  Zk = 2.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(3) 
  Zk = 6.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(4) 
  Zk = 8.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(5) 
  Zk = 10.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(6) 
  Zk = 11.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(7) 
  Zk = 12.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(8) 
  Zk = 14.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(9) 
  Zk = 15.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(10) 
  Zk = 16.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(11) 
  Zk = 18.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 
  Yk = Y(12) 
  Zk = 20.0_rp 
  ye = ye + Yk*Zk 
  iabar = iabar + Yk 
  fw = fw + Zk*Zk*Yk 

  abar = rp1/iabar 
  zbar = ye*abar 
  fw = fw + ye 

  iT = rp1/T 
  rhoye = rho*ye 
  gammap = CONST_S1*rhoye**othird*iT 
  gammapo4 = gammap**oquart 
  lngammap = log(gammap) 
  iT3 = iT*iT*iT 
  fw = CONST_S2*sqrt(fw*rho*iT3) 

  Z1 = 2.0_rp 
  Z2 = 6.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(11) = rates(11)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 8.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(12) = rates(12)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 10.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(13) = rates(13)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 11.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(14) = rates(14)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 12.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 24.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 24.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(15) = rates(15)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 14.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(16) = rates(16)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 15.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(17) = rates(17)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 16.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 32.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 32.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(18) = rates(18)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 18.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 36.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 36.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(19) = rates(19)*exp(H) 

  Z1 = 6.0_rp 
  Z2 = 6.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 12.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 12.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(20) = rates(20)*exp(H) 

  Z1 = 6.0_rp 
  Z2 = 6.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 12.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 12.0_rp 
   A2 = 12.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(21) = rates(21)*exp(H) 

  Z1 = 6.0_rp 
  Z2 = 8.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 12.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 12.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(22) = rates(22)*exp(H) 

  Z1 = 8.0_rp 
  Z2 = 8.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 16.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 16.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(23) = rates(23)*exp(H) 

  Z1 = 8.0_rp 
  Z2 = 8.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 16.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 16.0_rp 
   A2 = 16.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(24) = rates(24)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 10.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(25) = rates(25)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 10.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(26) = rates(26)*exp(H) 

  Z1 = 6.0_rp 
  Z2 = 10.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 12.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 12.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(27) = rates(27)*exp(H) 

  Z1 = 6.0_rp 
  Z2 = 10.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 12.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 12.0_rp 
   A2 = 20.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(28) = rates(28)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 11.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(29) = rates(29)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 11.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 23.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(30) = rates(30)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 12.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 24.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 24.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(31) = rates(31)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 14.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(32) = rates(32)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 14.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(33) = rates(33)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 14.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 28.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(34) = rates(34)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 15.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(35) = rates(35)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 15.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(36) = rates(36)*exp(H) 

  Z1 = 1.0_rp 
  Z2 = 15.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 1.0_rp 
   A2 = 31.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(37) = rates(37)*exp(H) 

  Z1 = 2.0_rp 
  Z2 = 2.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 4.0_rp 
   A2 = 4.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 4.0_rp 
   A2 = 4.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(38) = rates(38)*exp(H) 

  Z1 = 4.0_rp 
  Z2 = 2.0_rp 
  Zsum = Z1+Z2 
  Zprod = Z1*Z2 
  gammaeff = (rp2/Zsum)**othird*Zprod*gammap 
  if(gammaeff<0.3_rp) then 
   H = Zprod*fw 
  else if(gammaeff>0.8_rp) then 
   A1 = 8.0_rp 
   A2 = 4.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   H = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
  else 
   Hw = Zprod*fw 
   A1 = 8.0_rp 
   A2 = 4.0_rp 
   Ahat = A1*A2/(A1+A2) 
   Z1_5_12 = Z1**fvtwelfth 
   Z2_5_12 = Z2**fvtwelfth 
   Zsum_5_12 = Zsum**fvtwelfth 
   Z1_5_3 = Z1_5_12*Z1_5_12*Z1_5_12*Z1_5_12 
   Z2_5_3 = Z2_5_12*Z2_5_12*Z2_5_12*Z2_5_12 
   Zsum_5_3 = Zsum_5_12*Zsum_5_12*Zsum_5_12*Zsum_5_12 
   zt = Zsum_5_3 - Z1_5_3 - Z2_5_3 
   jt = Zsum_5_12 - Z1_5_12 - Z2_5_12 
   C = & 
   0.896434_rp*gammap*zt - & 
   3.44740_rp*gammapo4*jt - & 
   2.996_rp - & 
   0.5551_rp*(lngammap + fvthirds*log(Zprod/Zsum)) 
   tau = CONST_S3*(Zprod*Zprod*Ahat*iT)**othird 
   b = rp3*gammaeff/tau 
   b3 = b*b*b 
   b4 = b3*b 
   b5 = b4*b 
   b6 = b5*b 
   Hs = C - othird*tau*(rp5/rp32*b3-0.014_rp*b4-0.0128_rp*b5) - & 
   gammaeff*(0.0055_rp*b4-0.0098_rp*b5+0.0048_rp*b6) 
   H = rp2*(0.8_rp-gammaeff)*Hw + rp2*(gammaeff-0.3_rp)*Hs 
  end if 
  rates(38) = rates(38)*exp(H) 

end subroutine screen_rates 
#endif 

subroutine compute_network_residuals(Y,rho,R,res,jac,dedt,return_jac,return_dedt) 
  use source 
  real(kind=rp), dimension(1:nspecies), intent(in) :: Y 
  real(kind=rp), intent(in) :: rho 
  real(kind=rp), dimension(1:nreacs), intent(in) :: R 
  real(kind=rp), dimension(1:nspecies), intent(inout) :: res 
  real(kind=rp), dimension(1:nspecies,1:nspecies), intent(inout) :: jac 
  real(kind=rp), dimension(1:nreacs), intent(inout) :: dedt 
  logical, intent(in) :: return_jac 
  logical, intent(in) :: return_dedt 

  real(kind=rp) :: & 
  ye, & 
  lam1, & 
  lam2, & 
  lam3, & 
  lam4, & 
  lam5, & 
  lam6, & 
  lam7, & 
  lam8, & 
  lam9, & 
  lam10, & 
  lam11, & 
  lam12, & 
  lam13, & 
  lam14, & 
  lam15, & 
  lam16, & 
  lam17, & 
  lam18, & 
  lam19, & 
  lam20, & 
  lam21, & 
  lam22, & 
  lam23, & 
  lam24, & 
  lam25, & 
  lam26, & 
  lam27, & 
  lam28, & 
  lam29, & 
  lam30, & 
  lam31, & 
  lam32, & 
  lam33, & 
  lam34, & 
  lam35, & 
  lam36, & 
  lam37, & 
  lam38, & 
  R1, & 
  R2, & 
  R3, & 
  R4, & 
  R5, & 
  R6, & 
  R7, & 
  R8, & 
  R9, & 
  R10, & 
  R11, & 
  R12, & 
  R13, & 
  R14, & 
  R15, & 
  R16, & 
  R17, & 
  R18, & 
  R19, & 
  R20, & 
  R21, & 
  R22, & 
  R23, & 
  R24, & 
  R25, & 
  R26, & 
  R27, & 
  R28, & 
  R29, & 
  R30, & 
  R31, & 
  R32, & 
  R33, & 
  R34, & 
  R35, & 
  R36, & 
  R37, & 
  R38, & 
  Y1, & 
  Y2, & 
  Y3, & 
  Y4, & 
  Y5, & 
  Y6, & 
  Y7, & 
  Y8, & 
  Y9, & 
  Y10, & 
  Y11, & 
  Y12, & 
  dlam_dY_1_4, & 
  dlam_dY_2_5, & 
  dlam_dY_3_7, & 
  dlam_dY_4_7, & 
  dlam_dY_5_8, & 
  dlam_dY_6_10, & 
  dlam_dY_7_10, & 
  dlam_dY_8_11, & 
  dlam_dY_9_12, & 
  dlam_dY_10_3, & 
  dlam_dY_11_2, & 
  dlam_dY_11_3, & 
  dlam_dY_12_2, & 
  dlam_dY_12_4, & 
  dlam_dY_13_2, & 
  dlam_dY_13_5, & 
  dlam_dY_14_1, & 
  dlam_dY_14_6, & 
  dlam_dY_15_2, & 
  dlam_dY_15_7, & 
  dlam_dY_16_2, & 
  dlam_dY_16_8, & 
  dlam_dY_17_1, & 
  dlam_dY_17_9, & 
  dlam_dY_18_2, & 
  dlam_dY_18_10, & 
  dlam_dY_19_2, & 
  dlam_dY_19_11, & 
  dlam_dY_20_3, & 
  dlam_dY_21_3, & 
  dlam_dY_22_3, & 
  dlam_dY_22_4, & 
  dlam_dY_23_4, & 
  dlam_dY_24_4, & 
  dlam_dY_25_2, & 
  dlam_dY_25_5, & 
  dlam_dY_26_2, & 
  dlam_dY_26_5, & 
  dlam_dY_27_3, & 
  dlam_dY_27_5, & 
  dlam_dY_28_3, & 
  dlam_dY_28_5, & 
  dlam_dY_29_1, & 
  dlam_dY_29_6, & 
  dlam_dY_30_1, & 
  dlam_dY_30_6, & 
  dlam_dY_31_2, & 
  dlam_dY_31_7, & 
  dlam_dY_32_2, & 
  dlam_dY_32_8, & 
  dlam_dY_33_2, & 
  dlam_dY_33_8, & 
  dlam_dY_34_2, & 
  dlam_dY_34_8, & 
  dlam_dY_35_1, & 
  dlam_dY_35_9, & 
  dlam_dY_36_1, & 
  dlam_dY_36_9, & 
  dlam_dY_37_1, & 
  dlam_dY_37_9, & 
  dlam_dY_38_2, & 
  rho2 

  rho2 = rho*rho 
  ye=rp1 

  R1=R(1) 
  R2=R(2) 
  R3=R(3) 
  R4=R(4) 
  R5=R(5) 
  R6=R(6) 
  R7=R(7) 
  R8=R(8) 
  R9=R(9) 
  R10=R(10) 
  R11=R(11) 
  R12=R(12) 
  R13=R(13) 
  R14=R(14) 
  R15=R(15) 
  R16=R(16) 
  R17=R(17) 
  R18=R(18) 
  R19=R(19) 
  R20=R(20) 
  R21=R(21) 
  R22=R(22) 
  R23=R(23) 
  R24=R(24) 
  R25=R(25) 
  R26=R(26) 
  R27=R(27) 
  R28=R(28) 
  R29=R(29) 
  R30=R(30) 
  R31=R(31) 
  R32=R(32) 
  R33=R(33) 
  R34=R(34) 
  R35=R(35) 
  R36=R(36) 
  R37=R(37) 
  R38=R(38) 

  Y1=Y(1) 
  Y2=Y(2) 
  Y3=Y(3) 
  Y4=Y(4) 
  Y5=Y(5) 
  Y6=Y(6) 
  Y7=Y(7) 
  Y8=Y(8) 
  Y9=Y(9) 
  Y10=Y(10) 
  Y11=Y(11) 
  Y12=Y(12) 

  lam1=Y4*R1
  lam2=Y5*R2
  lam3=Y7*R3
  lam4=Y7*R4
  lam5=Y8*R5
  lam6=Y10*R6
  lam7=Y10*R7
  lam8=Y11*R8
  lam9=Y12*R9
  lam10=Y3*R10
  lam11=rho*Y2*Y3*R11
  lam12=rho*Y2*Y4*R12
  lam13=rho*Y2*Y5*R13
  lam14=rho*Y1*Y6*R14
  lam15=rho*Y2*Y7*R15
  lam16=rho*Y2*Y8*R16
  lam17=rho*Y1*Y9*R17
  lam18=rho*Y2*Y10*R18
  lam19=rho*Y2*Y11*R19
  lam20=rph*rho*Y3*Y3*R20
  lam21=rph*rho*Y3*Y3*R21
  lam22=rho*Y3*Y4*R22
  lam23=rph*rho*Y4*Y4*R23
  lam24=rph*rho*Y4*Y4*R24
  lam25=rho*Y2*Y5*R25
  lam26=rho*Y2*Y5*R26
  lam27=rho*Y3*Y5*R27
  lam28=rho*Y3*Y5*R28
  lam29=rho*Y1*Y6*R29
  lam30=rho*Y1*Y6*R30
  lam31=rho*Y2*Y7*R31
  lam32=rho*Y2*Y8*R32
  lam33=rho*Y2*Y8*R33
  lam34=rho*Y2*Y8*R34
  lam35=rho*Y1*Y9*R35
  lam36=rho*Y1*Y9*R36
  lam37=rho*Y1*Y9*R37
  lam38=osixth*rho2*Y2*Y2*Y2*R38

  res(1)= & 
  -lam3 & 
  -lam6 & 
  +lam14 & 
  +lam17 & 
  -lam20 & 
  -lam23 & 
  -lam25 & 
  -lam27 & 
  +lam29 & 
  +lam30 & 
  -lam32 & 
  +lam35 & 
  +lam36 & 
  +lam37

  res(2)= & 
  -lam1 & 
  -lam2 & 
  -lam4 & 
  -lam5 & 
  -lam7 & 
  -lam8 & 
  -lam9 & 
  -rp3*lam10 & 
  +lam11 & 
  +lam12 & 
  +lam13 & 
  +lam15 & 
  +lam16 & 
  +lam18 & 
  +lam19 & 
  -lam21 & 
  -lam22 & 
  -lam24 & 
  +lam25 & 
  +lam26 & 
  -lam28 & 
  -lam29 & 
  +lam31 & 
  +lam32 & 
  +lam33 & 
  +lam34 & 
  -lam35 & 
  +rp3*lam38

  res(3)= & 
  -lam1 & 
  +lam10 & 
  +lam11 & 
  +rp2*lam20 & 
  +rp2*lam21 & 
  +lam22 & 
  -rp2*lam26 & 
  +lam27 & 
  +lam28 & 
  -rp2*lam30 & 
  -lam31 & 
  -lam33 & 
  -lam36 & 
  -lam38

  res(4)= & 
  +lam1 & 
  -lam2 & 
  -lam11 & 
  +lam12 & 
  +lam22 & 
  +rp2*lam23 & 
  +rp2*lam24 & 
  -lam31 & 
  -rp2*lam34 & 
  -rp2*lam37

  res(5)= & 
  +lam2 & 
  -lam4 & 
  -lam12 & 
  +lam13 & 
  -lam21 & 
  +lam25 & 
  +lam26 & 
  +lam27 & 
  +lam28 & 
  -lam29 & 
  -lam33 & 
  -lam36

  res(6)= & 
  -lam3 & 
  +lam14 & 
  -lam20 & 
  -lam25 & 
  +lam29 & 
  +lam30

  res(7)= & 
  +lam3 & 
  +lam4 & 
  -lam5 & 
  -lam13 & 
  -lam14 & 
  +lam15 & 
  -lam22 & 
  +lam31

  res(8)= & 
  +lam5 & 
  -lam7 & 
  -lam15 & 
  +lam16 & 
  -lam24 & 
  -lam28 & 
  +lam32 & 
  +lam33 & 
  +lam34 & 
  -lam35

  res(9)= & 
  -lam6 & 
  +lam17 & 
  -lam23 & 
  -lam27 & 
  -lam32 & 
  +lam35 & 
  +lam36 & 
  +lam37

  res(10)= & 
  +lam6 & 
  +lam7 & 
  -lam8 & 
  -lam16 & 
  -lam17 & 
  +lam18

  res(11)= & 
  +lam8 & 
  -lam9 & 
  -lam18 & 
  +lam19

  res(12)= & 
  +lam9 & 
  -lam19

  if(return_jac) then 

   dlam_dY_1_4=R1 
   dlam_dY_2_5=R2 
   dlam_dY_3_7=R3 
   dlam_dY_4_7=R4 
   dlam_dY_5_8=R5 
   dlam_dY_6_10=R6 
   dlam_dY_7_10=R7 
   dlam_dY_8_11=R8 
   dlam_dY_9_12=R9 
   dlam_dY_10_3=R10 
   dlam_dY_11_2=rho*Y3*R11 
   dlam_dY_11_3=rho*Y2*R11 
   dlam_dY_12_2=rho*Y4*R12 
   dlam_dY_12_4=rho*Y2*R12 
   dlam_dY_13_2=rho*Y5*R13 
   dlam_dY_13_5=rho*Y2*R13 
   dlam_dY_14_1=rho*Y6*R14 
   dlam_dY_14_6=rho*Y1*R14 
   dlam_dY_15_2=rho*Y7*R15 
   dlam_dY_15_7=rho*Y2*R15 
   dlam_dY_16_2=rho*Y8*R16 
   dlam_dY_16_8=rho*Y2*R16 
   dlam_dY_17_1=rho*Y9*R17 
   dlam_dY_17_9=rho*Y1*R17 
   dlam_dY_18_2=rho*Y10*R18 
   dlam_dY_18_10=rho*Y2*R18 
   dlam_dY_19_2=rho*Y11*R19 
   dlam_dY_19_11=rho*Y2*R19 
   dlam_dY_20_3=rho*Y3*R20 
   dlam_dY_21_3=rho*Y3*R21 
   dlam_dY_22_3=rho*Y4*R22 
   dlam_dY_22_4=rho*Y3*R22 
   dlam_dY_23_4=rho*Y4*R23 
   dlam_dY_24_4=rho*Y4*R24 
   dlam_dY_25_2=rho*Y5*R25 
   dlam_dY_25_5=rho*Y2*R25 
   dlam_dY_26_2=rho*Y5*R26 
   dlam_dY_26_5=rho*Y2*R26 
   dlam_dY_27_3=rho*Y5*R27 
   dlam_dY_27_5=rho*Y3*R27 
   dlam_dY_28_3=rho*Y5*R28 
   dlam_dY_28_5=rho*Y3*R28 
   dlam_dY_29_1=rho*Y6*R29 
   dlam_dY_29_6=rho*Y1*R29 
   dlam_dY_30_1=rho*Y6*R30 
   dlam_dY_30_6=rho*Y1*R30 
   dlam_dY_31_2=rho*Y7*R31 
   dlam_dY_31_7=rho*Y2*R31 
   dlam_dY_32_2=rho*Y8*R32 
   dlam_dY_32_8=rho*Y2*R32 
   dlam_dY_33_2=rho*Y8*R33 
   dlam_dY_33_8=rho*Y2*R33 
   dlam_dY_34_2=rho*Y8*R34 
   dlam_dY_34_8=rho*Y2*R34 
   dlam_dY_35_1=rho*Y9*R35 
   dlam_dY_35_9=rho*Y1*R35 
   dlam_dY_36_1=rho*Y9*R36 
   dlam_dY_36_9=rho*Y1*R36 
   dlam_dY_37_1=rho*Y9*R37 
   dlam_dY_37_9=rho*Y1*R37 
   dlam_dY_38_2=rph*rho2*Y2*Y2*R38 

   jac(1,1)= & 
   +dlam_dY_14_1 & 
   +dlam_dY_17_1 & 
   +dlam_dY_29_1 & 
   +dlam_dY_30_1 & 
   +dlam_dY_35_1 & 
   +dlam_dY_36_1 & 
   +dlam_dY_37_1

   jac(1,2)= & 
   -dlam_dY_25_2 & 
   -dlam_dY_32_2

   jac(1,3)= & 
   -dlam_dY_20_3 & 
   -dlam_dY_27_3

   jac(1,4)= & 
   -dlam_dY_23_4

   jac(1,5)= & 
   -dlam_dY_25_5 & 
   -dlam_dY_27_5

   jac(1,6)= & 
   +dlam_dY_14_6 & 
   +dlam_dY_29_6 & 
   +dlam_dY_30_6

   jac(1,7)= & 
   -dlam_dY_3_7

   jac(1,8)= & 
   -dlam_dY_32_8

   jac(1,9)= & 
   +dlam_dY_17_9 & 
   +dlam_dY_35_9 & 
   +dlam_dY_36_9 & 
   +dlam_dY_37_9

   jac(1,10)= & 
   -dlam_dY_6_10

   jac(2,1)= & 
   -dlam_dY_29_1 & 
   -dlam_dY_35_1

   jac(2,2)= & 
   +dlam_dY_11_2 & 
   +dlam_dY_12_2 & 
   +dlam_dY_13_2 & 
   +dlam_dY_15_2 & 
   +dlam_dY_16_2 & 
   +dlam_dY_18_2 & 
   +dlam_dY_19_2 & 
   +dlam_dY_25_2 & 
   +dlam_dY_26_2 & 
   +dlam_dY_31_2 & 
   +dlam_dY_32_2 & 
   +dlam_dY_33_2 & 
   +dlam_dY_34_2 & 
   +rp3*dlam_dY_38_2

   jac(2,3)= & 
   -rp3*dlam_dY_10_3 & 
   +dlam_dY_11_3 & 
   -dlam_dY_21_3 & 
   -dlam_dY_22_3 & 
   -dlam_dY_28_3

   jac(2,4)= & 
   -dlam_dY_1_4 & 
   +dlam_dY_12_4 & 
   -dlam_dY_22_4 & 
   -dlam_dY_24_4

   jac(2,5)= & 
   -dlam_dY_2_5 & 
   +dlam_dY_13_5 & 
   +dlam_dY_25_5 & 
   +dlam_dY_26_5 & 
   -dlam_dY_28_5

   jac(2,6)= & 
   -dlam_dY_29_6

   jac(2,7)= & 
   -dlam_dY_4_7 & 
   +dlam_dY_15_7 & 
   +dlam_dY_31_7

   jac(2,8)= & 
   -dlam_dY_5_8 & 
   +dlam_dY_16_8 & 
   +dlam_dY_32_8 & 
   +dlam_dY_33_8 & 
   +dlam_dY_34_8

   jac(2,9)= & 
   -dlam_dY_35_9

   jac(2,10)= & 
   -dlam_dY_7_10 & 
   +dlam_dY_18_10

   jac(2,11)= & 
   -dlam_dY_8_11 & 
   +dlam_dY_19_11

   jac(2,12)= & 
   -dlam_dY_9_12

   jac(3,1)= & 
   -rp2*dlam_dY_30_1 & 
   -dlam_dY_36_1

   jac(3,2)= & 
   +dlam_dY_11_2 & 
   -rp2*dlam_dY_26_2 & 
   -dlam_dY_31_2 & 
   -dlam_dY_33_2 & 
   -dlam_dY_38_2

   jac(3,3)= & 
   +dlam_dY_10_3 & 
   +dlam_dY_11_3 & 
   +rp2*dlam_dY_20_3 & 
   +rp2*dlam_dY_21_3 & 
   +dlam_dY_22_3 & 
   +dlam_dY_27_3 & 
   +dlam_dY_28_3

   jac(3,4)= & 
   -dlam_dY_1_4 & 
   +dlam_dY_22_4

   jac(3,5)= & 
   -rp2*dlam_dY_26_5 & 
   +dlam_dY_27_5 & 
   +dlam_dY_28_5

   jac(3,6)= & 
   -rp2*dlam_dY_30_6

   jac(3,7)= & 
   -dlam_dY_31_7

   jac(3,8)= & 
   -dlam_dY_33_8

   jac(3,9)= & 
   -dlam_dY_36_9

   jac(4,1)= & 
   -rp2*dlam_dY_37_1

   jac(4,2)= & 
   -dlam_dY_11_2 & 
   +dlam_dY_12_2 & 
   -dlam_dY_31_2 & 
   -rp2*dlam_dY_34_2

   jac(4,3)= & 
   -dlam_dY_11_3 & 
   +dlam_dY_22_3

   jac(4,4)= & 
   +dlam_dY_1_4 & 
   +dlam_dY_12_4 & 
   +dlam_dY_22_4 & 
   +rp2*dlam_dY_23_4 & 
   +rp2*dlam_dY_24_4

   jac(4,5)= & 
   -dlam_dY_2_5

   jac(4,7)= & 
   -dlam_dY_31_7

   jac(4,8)= & 
   -rp2*dlam_dY_34_8

   jac(4,9)= & 
   -rp2*dlam_dY_37_9

   jac(5,1)= & 
   -dlam_dY_29_1 & 
   -dlam_dY_36_1

   jac(5,2)= & 
   -dlam_dY_12_2 & 
   +dlam_dY_13_2 & 
   +dlam_dY_25_2 & 
   +dlam_dY_26_2 & 
   -dlam_dY_33_2

   jac(5,3)= & 
   -dlam_dY_21_3 & 
   +dlam_dY_27_3 & 
   +dlam_dY_28_3

   jac(5,4)= & 
   -dlam_dY_12_4

   jac(5,5)= & 
   +dlam_dY_2_5 & 
   +dlam_dY_13_5 & 
   +dlam_dY_25_5 & 
   +dlam_dY_26_5 & 
   +dlam_dY_27_5 & 
   +dlam_dY_28_5

   jac(5,6)= & 
   -dlam_dY_29_6

   jac(5,7)= & 
   -dlam_dY_4_7

   jac(5,8)= & 
   -dlam_dY_33_8

   jac(5,9)= & 
   -dlam_dY_36_9

   jac(6,1)= & 
   +dlam_dY_14_1 & 
   +dlam_dY_29_1 & 
   +dlam_dY_30_1

   jac(6,2)= & 
   -dlam_dY_25_2

   jac(6,3)= & 
   -dlam_dY_20_3

   jac(6,5)= & 
   -dlam_dY_25_5

   jac(6,6)= & 
   +dlam_dY_14_6 & 
   +dlam_dY_29_6 & 
   +dlam_dY_30_6

   jac(6,7)= & 
   -dlam_dY_3_7

   jac(7,1)= & 
   -dlam_dY_14_1

   jac(7,2)= & 
   -dlam_dY_13_2 & 
   +dlam_dY_15_2 & 
   +dlam_dY_31_2

   jac(7,3)= & 
   -dlam_dY_22_3

   jac(7,4)= & 
   -dlam_dY_22_4

   jac(7,5)= & 
   -dlam_dY_13_5

   jac(7,6)= & 
   -dlam_dY_14_6

   jac(7,7)= & 
   +dlam_dY_3_7 & 
   +dlam_dY_4_7 & 
   +dlam_dY_15_7 & 
   +dlam_dY_31_7

   jac(7,8)= & 
   -dlam_dY_5_8

   jac(8,1)= & 
   -dlam_dY_35_1

   jac(8,2)= & 
   -dlam_dY_15_2 & 
   +dlam_dY_16_2 & 
   +dlam_dY_32_2 & 
   +dlam_dY_33_2 & 
   +dlam_dY_34_2

   jac(8,3)= & 
   -dlam_dY_28_3

   jac(8,4)= & 
   -dlam_dY_24_4

   jac(8,5)= & 
   -dlam_dY_28_5

   jac(8,7)= & 
   -dlam_dY_15_7

   jac(8,8)= & 
   +dlam_dY_5_8 & 
   +dlam_dY_16_8 & 
   +dlam_dY_32_8 & 
   +dlam_dY_33_8 & 
   +dlam_dY_34_8

   jac(8,9)= & 
   -dlam_dY_35_9

   jac(8,10)= & 
   -dlam_dY_7_10

   jac(9,1)= & 
   +dlam_dY_17_1 & 
   +dlam_dY_35_1 & 
   +dlam_dY_36_1 & 
   +dlam_dY_37_1

   jac(9,2)= & 
   -dlam_dY_32_2

   jac(9,3)= & 
   -dlam_dY_27_3

   jac(9,4)= & 
   -dlam_dY_23_4

   jac(9,5)= & 
   -dlam_dY_27_5

   jac(9,8)= & 
   -dlam_dY_32_8

   jac(9,9)= & 
   +dlam_dY_17_9 & 
   +dlam_dY_35_9 & 
   +dlam_dY_36_9 & 
   +dlam_dY_37_9

   jac(9,10)= & 
   -dlam_dY_6_10

   jac(10,1)= & 
   -dlam_dY_17_1

   jac(10,2)= & 
   -dlam_dY_16_2 & 
   +dlam_dY_18_2

   jac(10,8)= & 
   -dlam_dY_16_8

   jac(10,9)= & 
   -dlam_dY_17_9

   jac(10,10)= & 
   +dlam_dY_6_10 & 
   +dlam_dY_7_10 & 
   +dlam_dY_18_10

   jac(10,11)= & 
   -dlam_dY_8_11

   jac(11,2)= & 
   -dlam_dY_18_2 & 
   +dlam_dY_19_2

   jac(11,10)= & 
   -dlam_dY_18_10

   jac(11,11)= & 
   +dlam_dY_8_11 & 
   +dlam_dY_19_11

   jac(11,12)= & 
   -dlam_dY_9_12

   jac(12,2)= & 
   -dlam_dY_19_2

   jac(12,11)= & 
   -dlam_dY_19_11

   jac(12,12)= & 
   +dlam_dY_9_12

  end if 

  if(return_dedt) then 

   dedt(1)=-6.910202e+18_rp*lam1*rho 
   dedt(2)=-4.563611e+18_rp*lam2*rho 
   dedt(3)=-1.128174e+19_rp*lam3*rho 
   dedt(4)=-8.989104e+18_rp*lam4*rho 
   dedt(5)=-9.633230e+18_rp*lam5*rho 
   dedt(6)=-8.552247e+18_rp*lam6*rho 
   dedt(7)=-6.703627e+18_rp*lam7*rho 
   dedt(8)=-6.407359e+18_rp*lam8*rho 
   dedt(9)=-6.792230e+18_rp*lam9*rho 
   dedt(10)=-7.019308e+18_rp*lam10*rho 
   dedt(11)=6.910202e+18_rp*lam11*rho 
   dedt(12)=4.563611e+18_rp*lam12*rho 
   dedt(13)=8.989104e+18_rp*lam13*rho 
   dedt(14)=1.128174e+19_rp*lam14*rho 
   dedt(15)=9.633230e+18_rp*lam15*rho 
   dedt(16)=6.703627e+18_rp*lam16*rho 
   dedt(17)=8.552247e+18_rp*lam17*rho 
   dedt(18)=6.407359e+18_rp*lam18*rho 
   dedt(19)=6.792230e+18_rp*lam19*rho 
   dedt(20)=2.163201e+18_rp*lam20*rho 
   dedt(21)=4.458587e+18_rp*lam21*rho 
   dedt(22)=6.533022e+18_rp*lam22*rho 
   dedt(23)=7.408144e+18_rp*lam23*rho 
   dedt(24)=9.255838e+18_rp*lam24*rho 
   dedt(25)=-2.292617e+18_rp*lam25*rho 
   dedt(26)=-4.458587e+18_rp*lam26*rho 
   dedt(27)=9.748878e+18_rp*lam27*rho 
   dedt(28)=1.159860e+19_rp*lam28*rho 
   dedt(29)=2.292617e+18_rp*lam29*rho 
   dedt(30)=-2.163201e+18_rp*lam30*rho 
   dedt(31)=-6.533022e+18_rp*lam31*rho 
   dedt(32)=-1.848630e+18_rp*lam32*rho 
   dedt(33)=-1.159860e+19_rp*lam33*rho 
   dedt(34)=-9.255838e+18_rp*lam34*rho 
   dedt(35)=1.848630e+18_rp*lam35*rho 
   dedt(36)=-9.748878e+18_rp*lam36*rho 
   dedt(37)=-7.408144e+18_rp*lam37*rho 
   dedt(38)=7.019308e+18_rp*lam38*rho 

  end if 

end subroutine compute_network_residuals 
#endif 

#ifdef SAVE_SPECIES_FLUXES 
subroutine species_residuals_per_reac(Y,rho,R,Xds) 
  use source 
  real(kind=rp), dimension(1:nspecies), intent(in) :: Y 
  real(kind=rp), intent(in) :: rho 
  real(kind=rp), dimension(1:nreacs), intent(in) :: R 
  real(kind=rp), dimension(1:nspecies,1:nreacs), intent(inout) :: Xds 

  real(kind=rp) :: & 
  ye, & 
  lam1, & 
  lam2, & 
  lam3, & 
  lam4, & 
  lam5, & 
  lam6, & 
  lam7, & 
  lam8, & 
  lam9, & 
  lam10, & 
  lam11, & 
  lam12, & 
  lam13, & 
  lam14, & 
  lam15, & 
  lam16, & 
  lam17, & 
  lam18, & 
  lam19, & 
  lam20, & 
  lam21, & 
  lam22, & 
  lam23, & 
  lam24, & 
  lam25, & 
  lam26, & 
  lam27, & 
  lam28, & 
  lam29, & 
  lam30, & 
  lam31, & 
  lam32, & 
  lam33, & 
  lam34, & 
  lam35, & 
  lam36, & 
  lam37, & 
  lam38, & 
  R1, & 
  R2, & 
  R3, & 
  R4, & 
  R5, & 
  R6, & 
  R7, & 
  R8, & 
  R9, & 
  R10, & 
  R11, & 
  R12, & 
  R13, & 
  R14, & 
  R15, & 
  R16, & 
  R17, & 
  R18, & 
  R19, & 
  R20, & 
  R21, & 
  R22, & 
  R23, & 
  R24, & 
  R25, & 
  R26, & 
  R27, & 
  R28, & 
  R29, & 
  R30, & 
  R31, & 
  R32, & 
  R33, & 
  R34, & 
  R35, & 
  R36, & 
  R37, & 
  R38, & 
  Y1, & 
  Y2, & 
  Y3, & 
  Y4, & 
  Y5, & 
  Y6, & 
  Y7, & 
  Y8, & 
  Y9, & 
  Y10, & 
  Y11, & 
  Y12, & 
  rho2 

  rho2 = rho*rho 
  ye=rp1 

  R1=R(1) 
  R2=R(2) 
  R3=R(3) 
  R4=R(4) 
  R5=R(5) 
  R6=R(6) 
  R7=R(7) 
  R8=R(8) 
  R9=R(9) 
  R10=R(10) 
  R11=R(11) 
  R12=R(12) 
  R13=R(13) 
  R14=R(14) 
  R15=R(15) 
  R16=R(16) 
  R17=R(17) 
  R18=R(18) 
  R19=R(19) 
  R20=R(20) 
  R21=R(21) 
  R22=R(22) 
  R23=R(23) 
  R24=R(24) 
  R25=R(25) 
  R26=R(26) 
  R27=R(27) 
  R28=R(28) 
  R29=R(29) 
  R30=R(30) 
  R31=R(31) 
  R32=R(32) 
  R33=R(33) 
  R34=R(34) 
  R35=R(35) 
  R36=R(36) 
  R37=R(37) 
  R38=R(38) 

  Y1=Y(1) 
  Y2=Y(2) 
  Y3=Y(3) 
  Y4=Y(4) 
  Y5=Y(5) 
  Y6=Y(6) 
  Y7=Y(7) 
  Y8=Y(8) 
  Y9=Y(9) 
  Y10=Y(10) 
  Y11=Y(11) 
  Y12=Y(12) 

  lam1=Y4*R1
  lam2=Y5*R2
  lam3=Y7*R3
  lam4=Y7*R4
  lam5=Y8*R5
  lam6=Y10*R6
  lam7=Y10*R7
  lam8=Y11*R8
  lam9=Y12*R9
  lam10=Y3*R10
  lam11=rho*Y2*Y3*R11
  lam12=rho*Y2*Y4*R12
  lam13=rho*Y2*Y5*R13
  lam14=rho*Y1*Y6*R14
  lam15=rho*Y2*Y7*R15
  lam16=rho*Y2*Y8*R16
  lam17=rho*Y1*Y9*R17
  lam18=rho*Y2*Y10*R18
  lam19=rho*Y2*Y11*R19
  lam20=rph*rho*Y3*Y3*R20
  lam21=rph*rho*Y3*Y3*R21
  lam22=rho*Y3*Y4*R22
  lam23=rph*rho*Y4*Y4*R23
  lam24=rph*rho*Y4*Y4*R24
  lam25=rho*Y2*Y5*R25
  lam26=rho*Y2*Y5*R26
  lam27=rho*Y3*Y5*R27
  lam28=rho*Y3*Y5*R28
  lam29=rho*Y1*Y6*R29
  lam30=rho*Y1*Y6*R30
  lam31=rho*Y2*Y7*R31
  lam32=rho*Y2*Y8*R32
  lam33=rho*Y2*Y8*R33
  lam34=rho*Y2*Y8*R34
  lam35=rho*Y1*Y9*R35
  lam36=rho*Y1*Y9*R36
  lam37=rho*Y1*Y9*R37
  lam38=osixth*rho2*Y2*Y2*Y2*R38

  Xds(1,1) = rp0 
  Xds(1,2) = rp0 
  Xds(1,3) = lam3 
  Xds(1,4) = rp0 
  Xds(1,5) = rp0 
  Xds(1,6) = lam6 
  Xds(1,7) = rp0 
  Xds(1,8) = rp0 
  Xds(1,9) = rp0 
  Xds(1,10) = rp0 
  Xds(1,11) = rp0 
  Xds(1,12) = rp0 
  Xds(1,13) = rp0 
  Xds(1,14) = -lam14 
  Xds(1,15) = rp0 
  Xds(1,16) = rp0 
  Xds(1,17) = -lam17 
  Xds(1,18) = rp0 
  Xds(1,19) = rp0 
  Xds(1,20) = lam20 
  Xds(1,21) = rp0 
  Xds(1,22) = rp0 
  Xds(1,23) = lam23 
  Xds(1,24) = rp0 
  Xds(1,25) = lam25 
  Xds(1,26) = rp0 
  Xds(1,27) = lam27 
  Xds(1,28) = rp0 
  Xds(1,29) = -lam29 
  Xds(1,30) = -lam30 
  Xds(1,31) = rp0 
  Xds(1,32) = lam32 
  Xds(1,33) = rp0 
  Xds(1,34) = rp0 
  Xds(1,35) = -lam35 
  Xds(1,36) = -lam36 
  Xds(1,37) = -lam37 
  Xds(1,38) = rp0 
  Xds(2,1) = lam1 
  Xds(2,2) = lam2 
  Xds(2,3) = rp0 
  Xds(2,4) = lam4 
  Xds(2,5) = lam5 
  Xds(2,6) = rp0 
  Xds(2,7) = lam7 
  Xds(2,8) = lam8 
  Xds(2,9) = lam9 
  Xds(2,10) = rp3*lam10 
  Xds(2,11) = -lam11 
  Xds(2,12) = -lam12 
  Xds(2,13) = -lam13 
  Xds(2,14) = rp0 
  Xds(2,15) = -lam15 
  Xds(2,16) = -lam16 
  Xds(2,17) = rp0 
  Xds(2,18) = -lam18 
  Xds(2,19) = -lam19 
  Xds(2,20) = rp0 
  Xds(2,21) = lam21 
  Xds(2,22) = lam22 
  Xds(2,23) = rp0 
  Xds(2,24) = lam24 
  Xds(2,25) = -lam25 
  Xds(2,26) = -lam26 
  Xds(2,27) = rp0 
  Xds(2,28) = lam28 
  Xds(2,29) = lam29 
  Xds(2,30) = rp0 
  Xds(2,31) = -lam31 
  Xds(2,32) = -lam32 
  Xds(2,33) = -lam33 
  Xds(2,34) = -lam34 
  Xds(2,35) = lam35 
  Xds(2,36) = rp0 
  Xds(2,37) = rp0 
  Xds(2,38) = -rp3*lam38 
  Xds(3,1) = lam1 
  Xds(3,2) = rp0 
  Xds(3,3) = rp0 
  Xds(3,4) = rp0 
  Xds(3,5) = rp0 
  Xds(3,6) = rp0 
  Xds(3,7) = rp0 
  Xds(3,8) = rp0 
  Xds(3,9) = rp0 
  Xds(3,10) = -lam10 
  Xds(3,11) = -lam11 
  Xds(3,12) = rp0 
  Xds(3,13) = rp0 
  Xds(3,14) = rp0 
  Xds(3,15) = rp0 
  Xds(3,16) = rp0 
  Xds(3,17) = rp0 
  Xds(3,18) = rp0 
  Xds(3,19) = rp0 
  Xds(3,20) = -rp2*lam20 
  Xds(3,21) = -rp2*lam21 
  Xds(3,22) = -lam22 
  Xds(3,23) = rp0 
  Xds(3,24) = rp0 
  Xds(3,25) = rp0 
  Xds(3,26) = rp2*lam26 
  Xds(3,27) = -lam27 
  Xds(3,28) = -lam28 
  Xds(3,29) = rp0 
  Xds(3,30) = rp2*lam30 
  Xds(3,31) = lam31 
  Xds(3,32) = rp0 
  Xds(3,33) = lam33 
  Xds(3,34) = rp0 
  Xds(3,35) = rp0 
  Xds(3,36) = lam36 
  Xds(3,37) = rp0 
  Xds(3,38) = lam38 
  Xds(4,1) = -lam1 
  Xds(4,2) = lam2 
  Xds(4,3) = rp0 
  Xds(4,4) = rp0 
  Xds(4,5) = rp0 
  Xds(4,6) = rp0 
  Xds(4,7) = rp0 
  Xds(4,8) = rp0 
  Xds(4,9) = rp0 
  Xds(4,10) = rp0 
  Xds(4,11) = lam11 
  Xds(4,12) = -lam12 
  Xds(4,13) = rp0 
  Xds(4,14) = rp0 
  Xds(4,15) = rp0 
  Xds(4,16) = rp0 
  Xds(4,17) = rp0 
  Xds(4,18) = rp0 
  Xds(4,19) = rp0 
  Xds(4,20) = rp0 
  Xds(4,21) = rp0 
  Xds(4,22) = -lam22 
  Xds(4,23) = -rp2*lam23 
  Xds(4,24) = -rp2*lam24 
  Xds(4,25) = rp0 
  Xds(4,26) = rp0 
  Xds(4,27) = rp0 
  Xds(4,28) = rp0 
  Xds(4,29) = rp0 
  Xds(4,30) = rp0 
  Xds(4,31) = lam31 
  Xds(4,32) = rp0 
  Xds(4,33) = rp0 
  Xds(4,34) = rp2*lam34 
  Xds(4,35) = rp0 
  Xds(4,36) = rp0 
  Xds(4,37) = rp2*lam37 
  Xds(4,38) = rp0 
  Xds(5,1) = rp0 
  Xds(5,2) = -lam2 
  Xds(5,3) = rp0 
  Xds(5,4) = lam4 
  Xds(5,5) = rp0 
  Xds(5,6) = rp0 
  Xds(5,7) = rp0 
  Xds(5,8) = rp0 
  Xds(5,9) = rp0 
  Xds(5,10) = rp0 
  Xds(5,11) = rp0 
  Xds(5,12) = lam12 
  Xds(5,13) = -lam13 
  Xds(5,14) = rp0 
  Xds(5,15) = rp0 
  Xds(5,16) = rp0 
  Xds(5,17) = rp0 
  Xds(5,18) = rp0 
  Xds(5,19) = rp0 
  Xds(5,20) = rp0 
  Xds(5,21) = lam21 
  Xds(5,22) = rp0 
  Xds(5,23) = rp0 
  Xds(5,24) = rp0 
  Xds(5,25) = -lam25 
  Xds(5,26) = -lam26 
  Xds(5,27) = -lam27 
  Xds(5,28) = -lam28 
  Xds(5,29) = lam29 
  Xds(5,30) = rp0 
  Xds(5,31) = rp0 
  Xds(5,32) = rp0 
  Xds(5,33) = lam33 
  Xds(5,34) = rp0 
  Xds(5,35) = rp0 
  Xds(5,36) = lam36 
  Xds(5,37) = rp0 
  Xds(5,38) = rp0 
  Xds(6,1) = rp0 
  Xds(6,2) = rp0 
  Xds(6,3) = lam3 
  Xds(6,4) = rp0 
  Xds(6,5) = rp0 
  Xds(6,6) = rp0 
  Xds(6,7) = rp0 
  Xds(6,8) = rp0 
  Xds(6,9) = rp0 
  Xds(6,10) = rp0 
  Xds(6,11) = rp0 
  Xds(6,12) = rp0 
  Xds(6,13) = rp0 
  Xds(6,14) = -lam14 
  Xds(6,15) = rp0 
  Xds(6,16) = rp0 
  Xds(6,17) = rp0 
  Xds(6,18) = rp0 
  Xds(6,19) = rp0 
  Xds(6,20) = lam20 
  Xds(6,21) = rp0 
  Xds(6,22) = rp0 
  Xds(6,23) = rp0 
  Xds(6,24) = rp0 
  Xds(6,25) = lam25 
  Xds(6,26) = rp0 
  Xds(6,27) = rp0 
  Xds(6,28) = rp0 
  Xds(6,29) = -lam29 
  Xds(6,30) = -lam30 
  Xds(6,31) = rp0 
  Xds(6,32) = rp0 
  Xds(6,33) = rp0 
  Xds(6,34) = rp0 
  Xds(6,35) = rp0 
  Xds(6,36) = rp0 
  Xds(6,37) = rp0 
  Xds(6,38) = rp0 
  Xds(7,1) = rp0 
  Xds(7,2) = rp0 
  Xds(7,3) = -lam3 
  Xds(7,4) = -lam4 
  Xds(7,5) = lam5 
  Xds(7,6) = rp0 
  Xds(7,7) = rp0 
  Xds(7,8) = rp0 
  Xds(7,9) = rp0 
  Xds(7,10) = rp0 
  Xds(7,11) = rp0 
  Xds(7,12) = rp0 
  Xds(7,13) = lam13 
  Xds(7,14) = lam14 
  Xds(7,15) = -lam15 
  Xds(7,16) = rp0 
  Xds(7,17) = rp0 
  Xds(7,18) = rp0 
  Xds(7,19) = rp0 
  Xds(7,20) = rp0 
  Xds(7,21) = rp0 
  Xds(7,22) = lam22 
  Xds(7,23) = rp0 
  Xds(7,24) = rp0 
  Xds(7,25) = rp0 
  Xds(7,26) = rp0 
  Xds(7,27) = rp0 
  Xds(7,28) = rp0 
  Xds(7,29) = rp0 
  Xds(7,30) = rp0 
  Xds(7,31) = -lam31 
  Xds(7,32) = rp0 
  Xds(7,33) = rp0 
  Xds(7,34) = rp0 
  Xds(7,35) = rp0 
  Xds(7,36) = rp0 
  Xds(7,37) = rp0 
  Xds(7,38) = rp0 
  Xds(8,1) = rp0 
  Xds(8,2) = rp0 
  Xds(8,3) = rp0 
  Xds(8,4) = rp0 
  Xds(8,5) = -lam5 
  Xds(8,6) = rp0 
  Xds(8,7) = lam7 
  Xds(8,8) = rp0 
  Xds(8,9) = rp0 
  Xds(8,10) = rp0 
  Xds(8,11) = rp0 
  Xds(8,12) = rp0 
  Xds(8,13) = rp0 
  Xds(8,14) = rp0 
  Xds(8,15) = lam15 
  Xds(8,16) = -lam16 
  Xds(8,17) = rp0 
  Xds(8,18) = rp0 
  Xds(8,19) = rp0 
  Xds(8,20) = rp0 
  Xds(8,21) = rp0 
  Xds(8,22) = rp0 
  Xds(8,23) = rp0 
  Xds(8,24) = lam24 
  Xds(8,25) = rp0 
  Xds(8,26) = rp0 
  Xds(8,27) = rp0 
  Xds(8,28) = lam28 
  Xds(8,29) = rp0 
  Xds(8,30) = rp0 
  Xds(8,31) = rp0 
  Xds(8,32) = -lam32 
  Xds(8,33) = -lam33 
  Xds(8,34) = -lam34 
  Xds(8,35) = lam35 
  Xds(8,36) = rp0 
  Xds(8,37) = rp0 
  Xds(8,38) = rp0 
  Xds(9,1) = rp0 
  Xds(9,2) = rp0 
  Xds(9,3) = rp0 
  Xds(9,4) = rp0 
  Xds(9,5) = rp0 
  Xds(9,6) = lam6 
  Xds(9,7) = rp0 
  Xds(9,8) = rp0 
  Xds(9,9) = rp0 
  Xds(9,10) = rp0 
  Xds(9,11) = rp0 
  Xds(9,12) = rp0 
  Xds(9,13) = rp0 
  Xds(9,14) = rp0 
  Xds(9,15) = rp0 
  Xds(9,16) = rp0 
  Xds(9,17) = -lam17 
  Xds(9,18) = rp0 
  Xds(9,19) = rp0 
  Xds(9,20) = rp0 
  Xds(9,21) = rp0 
  Xds(9,22) = rp0 
  Xds(9,23) = lam23 
  Xds(9,24) = rp0 
  Xds(9,25) = rp0 
  Xds(9,26) = rp0 
  Xds(9,27) = lam27 
  Xds(9,28) = rp0 
  Xds(9,29) = rp0 
  Xds(9,30) = rp0 
  Xds(9,31) = rp0 
  Xds(9,32) = lam32 
  Xds(9,33) = rp0 
  Xds(9,34) = rp0 
  Xds(9,35) = -lam35 
  Xds(9,36) = -lam36 
  Xds(9,37) = -lam37 
  Xds(9,38) = rp0 
  Xds(10,1) = rp0 
  Xds(10,2) = rp0 
  Xds(10,3) = rp0 
  Xds(10,4) = rp0 
  Xds(10,5) = rp0 
  Xds(10,6) = -lam6 
  Xds(10,7) = -lam7 
  Xds(10,8) = lam8 
  Xds(10,9) = rp0 
  Xds(10,10) = rp0 
  Xds(10,11) = rp0 
  Xds(10,12) = rp0 
  Xds(10,13) = rp0 
  Xds(10,14) = rp0 
  Xds(10,15) = rp0 
  Xds(10,16) = lam16 
  Xds(10,17) = lam17 
  Xds(10,18) = -lam18 
  Xds(10,19) = rp0 
  Xds(10,20) = rp0 
  Xds(10,21) = rp0 
  Xds(10,22) = rp0 
  Xds(10,23) = rp0 
  Xds(10,24) = rp0 
  Xds(10,25) = rp0 
  Xds(10,26) = rp0 
  Xds(10,27) = rp0 
  Xds(10,28) = rp0 
  Xds(10,29) = rp0 
  Xds(10,30) = rp0 
  Xds(10,31) = rp0 
  Xds(10,32) = rp0 
  Xds(10,33) = rp0 
  Xds(10,34) = rp0 
  Xds(10,35) = rp0 
  Xds(10,36) = rp0 
  Xds(10,37) = rp0 
  Xds(10,38) = rp0 
  Xds(11,1) = rp0 
  Xds(11,2) = rp0 
  Xds(11,3) = rp0 
  Xds(11,4) = rp0 
  Xds(11,5) = rp0 
  Xds(11,6) = rp0 
  Xds(11,7) = rp0 
  Xds(11,8) = -lam8 
  Xds(11,9) = lam9 
  Xds(11,10) = rp0 
  Xds(11,11) = rp0 
  Xds(11,12) = rp0 
  Xds(11,13) = rp0 
  Xds(11,14) = rp0 
  Xds(11,15) = rp0 
  Xds(11,16) = rp0 
  Xds(11,17) = rp0 
  Xds(11,18) = lam18 
  Xds(11,19) = -lam19 
  Xds(11,20) = rp0 
  Xds(11,21) = rp0 
  Xds(11,22) = rp0 
  Xds(11,23) = rp0 
  Xds(11,24) = rp0 
  Xds(11,25) = rp0 
  Xds(11,26) = rp0 
  Xds(11,27) = rp0 
  Xds(11,28) = rp0 
  Xds(11,29) = rp0 
  Xds(11,30) = rp0 
  Xds(11,31) = rp0 
  Xds(11,32) = rp0 
  Xds(11,33) = rp0 
  Xds(11,34) = rp0 
  Xds(11,35) = rp0 
  Xds(11,36) = rp0 
  Xds(11,37) = rp0 
  Xds(11,38) = rp0 
  Xds(12,1) = rp0 
  Xds(12,2) = rp0 
  Xds(12,3) = rp0 
  Xds(12,4) = rp0 
  Xds(12,5) = rp0 
  Xds(12,6) = rp0 
  Xds(12,7) = rp0 
  Xds(12,8) = rp0 
  Xds(12,9) = -lam9 
  Xds(12,10) = rp0 
  Xds(12,11) = rp0 
  Xds(12,12) = rp0 
  Xds(12,13) = rp0 
  Xds(12,14) = rp0 
  Xds(12,15) = rp0 
  Xds(12,16) = rp0 
  Xds(12,17) = rp0 
  Xds(12,18) = rp0 
  Xds(12,19) = lam19 
  Xds(12,20) = rp0 
  Xds(12,21) = rp0 
  Xds(12,22) = rp0 
  Xds(12,23) = rp0 
  Xds(12,24) = rp0 
  Xds(12,25) = rp0 
  Xds(12,26) = rp0 
  Xds(12,27) = rp0 
  Xds(12,28) = rp0 
  Xds(12,29) = rp0 
  Xds(12,30) = rp0 
  Xds(12,31) = rp0 
  Xds(12,32) = rp0 
  Xds(12,33) = rp0 
  Xds(12,34) = rp0 
  Xds(12,35) = rp0 
  Xds(12,36) = rp0 
  Xds(12,37) = rp0 
  Xds(12,38) = rp0 

end subroutine species_residuals_per_reac 
#endif 

