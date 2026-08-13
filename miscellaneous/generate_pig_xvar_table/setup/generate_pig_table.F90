program generate_pig_table
  use source
  implicit none
 
  type(mpigrid) :: mgrid
  type(teos), dimension(1:int(Nrho/ddx1+1),1:int(NT/ddx2+1)) :: table 
  type(teos), dimension(1:Nrho+1,1:NT+1) :: table_glob
  integer :: ierr,iv
  real(kind=rp) :: Zm

  call init_mpi(mgrid)
  
  mgrid%X(:) = 0.0_rp
  mgrid%X(i_c12) = 0.00214939_rp
  mgrid%X(i_n14) = 0.0033643_rp
  mgrid%X(i_o16) = 0.00839068_rp
  mgrid%X(i_ne20) = 0.00621659_rp

  Zm = 0.0_rp
  do iv=i_c12,i_fe56
   Zm = Zm + mgrid%X(iv) 
  end do

  mgrid%X(i_h1) = ((1.0_rp-Zm)/real(X_nint,kind=rp))*real(sim_index,kind=rp)
  mgrid%X(i_he4) = 1.0_rp-Zm-mgrid%X(i_h1)
 
  call mpi_barrier(mgrid%comm_cart,ierr)

  call create_table_briquette(mgrid,table)

  call create_global_table(mgrid,table_glob)

  call mpi_barrier(mgrid%comm_cart,ierr)

  call mpi_finalize(ierr)

end program generate_pig_table
