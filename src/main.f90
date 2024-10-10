module global_variables
  implicit none
  !mathematical parameters
  real(8),parameter:: pi=4d0*atan(1d0)
  complex(8),parameter:: zi=(0d0,1d0)

  !physical systems
  complex(8),allocatable:: zpsi(:,:),zpsi_gs(:,:) !zpsi(2,nk)
  real(8),allocatable:: eps_bk(:,:)
  integer:: nk,nk1,nk2
  integer:: nk_s,nk_e,nk_ave,nk_remainder
  real(8):: a_vec(2,2),a_lattice, b_vec(2,2)
  real(8)::delta_vec(2,3)
  real(8),allocatable:: kx0(:),ky0(:),kxt(:),kyt(:)
  real(8):: t0_hop,eps_b,eps_n

  !time propagation
  integer:: nt,p
  real(8):: dt,Tprop

   !laser dields
  real(8),allocatable:: Act(:,:)
  real(8):: E0,omega0,Tpulse0

   include 'mpif.h'
   integer:: ierr,Nprocs, Myrank
 

end module global_variables
!--------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
program main
  use global_variables
  implicit none
  call MPI_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Myrank,ierr)
  call input_variables
  call initialize
  call calc_ground_state
  !write(*,*)"ok00"
 ! call band_1d
call time_propagation
 
call MPI_FINALIZE(ierr)
 
end program main
!-----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------

subroutine read_input_variables(filename)
  use global_variables
  implicit none
  character(*), intent(in) :: filename
  integer :: ios
  character(100) :: line
  character(100) :: variable_name
  character(100) :: equals
  character(100) :: value_str ! To read the value as a string
  real(8) :: variable_value

  ! Open the input file
  open(unit=10, file=filename, status='old', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'Error: Unable to open input file ', filename
     stop
  end if

  ! Read parameters from the input file
  do
  read(10, '(A)', iostat=ios) line
  if (ios /= 0) exit

  ! Parse the line to get variable name, '=', and variable value
  read(line, *, iostat=ios) variable_name, equals, variable_value
  select case (trim(adjustl(variable_name)))
  case ("nk1")
      nk1 = variable_value
  case ("nk2")
      nk2 = variable_value
  case ("Tprop")
      Tprop = variable_value
  case ("dt")
      dt = variable_value
  case ("E0")
      E0 = variable_value
  case ("omega0")
      omega0 = variable_value
  case ("Tpulse0")
      Tpulse0 = variable_value
      ! Add cases for other variables...
  case default
      write(*, *) 'Warning: Unrecognized variable "', trim(adjustl(variable_name)), '"'
  end select
  end do

  ! Close the input file
  close(unit=10)
write(*, *) '=========================================='
write(*,*) 'nk1=', nk1
write(*,*) 'nk2=', nk2
write(*,*) 'Tprop=', Tprop
write(*,*) 'dt=', dt
write(*,*) 'E0=',E0 
write(*,*) 'omega0=', omega0
write(*,*) 'Tpulse0=', Tpulse0 
write(*, *) '=========================================='
end subroutine read_input_variables


subroutine input_variables
  use global_variables
  implicit none
  
  character(100) :: input_filename
  ! Initialize input_filename with the name of your input file
  input_filename = './inp'
  
  call read_input_variables(input_filename) 
  !physical parameters
  !eps_b=0d0
  !eps_n=0d0
  eps_b=3.34d0/27.2114d0
  eps_n=-2.56d0/27.2114d0
  t0_hop=2.64d0/27.2114d0
 ! t0_hop=1d0
 
  !lattice constant
  a_lattice=2.456d0/0.5291772d0  !!2.5AA
  
  !lattice vector
  a_vec(1,1)=a_lattice*sqrt(3d0)/2d0
  a_vec(2,1)=a_lattice*(0.5d0)

   a_vec(1,2)=a_lattice*sqrt(3d0)/2d0
  a_vec(2,2)=a_lattice*(-0.5d0)
  
  !desplacement
  delta_vec(1,2) = a_lattice/sqrt(3d0)
  delta_vec(2,2) = a_lattice*0d0

  delta_vec(1,1)=a_lattice*(-1d0/(2d0*sqrt(3d0)))
  delta_vec(2,1) = a_lattice*0.5d0

  delta_vec(1,3)=a_lattice*(-1d0/(2d0*sqrt(3d0)))
  delta_vec(2,3) = a_lattice*(-0.5d0)

  !time propagation
  nt=aint(Tprop/dt)+1
  
end subroutine input_variables
!-----------------------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  real(8)::volume
  integer::ik1,ik2,ik

  nk = nk1*nk2

  nk_ave=nk/Nprocs
  nk_remainder=mod(nk,Nprocs)

  if(myrank<nk_remainder)then
  nk_s=(nk_ave+1)*myrank+1
  nk_e=(nk_ave+1)*myrank +(nk_ave+1)
else
   nk_s=(myrank-nk_remainder)*nk_ave+(nk_ave+1)*nk_remainder+1
   nk_e=(myrank-nk_remainder)*nk_ave+(nk_ave+1)*nk_remainder+nk_ave
end if

if(myrank == 0)then
   write(*,*)"nk=",nk
 end if

   write(*,*) "nk_s,nk_e",nk_s,nk_e,myrank

 ! call MPI_FINALIZE(ierr)
 ! stop
  
  allocate(kx0(nk),ky0(nk),kxt(nk),kyt(nk))
  allocate(zpsi(2,nk), eps_bk(2,nk),zpsi_gs(2,nk))
  
   !reciprocal vectors
   volume=a_vec(1,1)*a_vec(2,2)-a_vec(2,1)*a_vec(1,2)

  b_vec(1,1)=2d0*pi*(a_vec(2,2)*1d0)/volume
  b_vec(2,1)=2d0*pi*(-a_vec(1,2)*1d0)/volume
  
  b_vec(1,2)= 2d0*pi*(-1d0*a_vec(2,1))/volume
  b_vec(2,2)= 2d0*pi*(1d0*a_vec(1,1))/volume
  
  !if (myrank==0)then
 ! write(*,*) b_vec
!write(*,*)sum(a_vec(:,1)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,1)*b_vec(:,2))/(2d0*pi)
!write(*,*)sum(a_vec(:,2)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,2)*b_vec(:,2))/(2d0*pi)  
!  end if
	!kx0(1)=b_vec(1,1)*2d0/3d0+b_vec(1,2)*1d0/3d0
        !ky0(1)=b_vec(2,1)*2d0/3d0+b_vec(2,2)*1d0/3d0
	!kx0(2)=b_vec(1,1)*1d0/3d0+b_vec(1,2)*2d0/3d0
        !ky0(2)=b_vec(2,1)*1d0/3d0+b_vec(2,2)*2d0/3d0	
  ik=0
 do ik1=0,nk1-1
   do ik2 =0,nk2-1
        ik=ik+1

      
       kx0(ik)=b_vec(1,1)*ik1/dble(nk1)+b_vec(1,2)*ik2/dble(nk2)
       ky0(ik)=b_vec(2,1)*ik1/dble(nk1)+b_vec(2,2)*ik2/dble(nk2)
        
     end do
  end do
end subroutine initialize
!--------------------------------------------------------------------------------------------------------------------------------
subroutine calc_ground_state
  use global_variables
  implicit none
integer:: ik,ik1,ik2
complex(8)::zham(2,2),zvec(2,2),zfk
real(8)::kx_t,ky_t
real(8):: eps_t(2)
real(8),allocatable:: eps_bk_l(:,:)

allocate(eps_bk_l(2,nk))
eps_bk_l=0d0
do ik=nk_s,nk_e
   kx_t=kx0(ik)
   ky_t=ky0(ik)

   zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))
   
   call calc_eig_vec_2x2(zham,zvec,eps_t)
   zpsi(:,ik)=zvec(:,1)
   eps_bk_l(:,ik)=eps_t
   zpsi_gs(:,ik)=zvec(:,1)
   
end do

call MPI_ALLREDUCE(eps_bk_l,eps_bk,2*nk,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myrank==0)then

!open(20,file='band_map.out')
ik=0
do ik1=0,nk1-1
 do ik2=0,nk2-1
   ik=ik+1
 !  write(20,"(999e26.16e3)")kx0(ik),ky0(ik),eps_bk(1:2,ik)
 end do
 !write(20,*)
end do
!close(20)

end if

!call MPI_FINALIZE(ierr)
!stop

end subroutine calc_ground_state

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine calc_eig_vec_2x2(zham, zvec, eps_t)
  implicit none
  complex(8),intent(in) :: zham(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out) :: eps_t(2)
  real(8) :: a,c
  complex(8):: zb, zx, zy

! (a    zb)
! (zb^*  c)

  a = zham(1,1)
  zb = zham(1,2)
  c = zham(2,2)

  eps_t(1) = 0.5d0*( (a+c)-sqrt((a-c)**2+4d0*abs(zb)**2))
  eps_t(2) = 0.5d0*( (a+c)+sqrt((a-c)**2+4d0*abs(zb)**2))

  if(a<c)then
    zy = conjg(zb)/(eps_t(1)-c)
    zvec(1,1) = 1d0/sqrt(1d0+abs(zy)**2)
    zvec(2,1) = zy/sqrt(1d0+abs(zy)**2)

    zx = zb/(eps_t(2)-a)
    zvec(1,2) = zx/sqrt(1d0+abs(zx)**2)
    zvec(2,2) = 1d0/sqrt(1d0+abs(zx)**2)

  else

    zx = zb/(eps_t(1)-a)
    zvec(1,1) = zx/sqrt(1d0+abs(zx)**2)
    zvec(2,1) = 1d0/sqrt(1d0+abs(zx)**2)

    zy = conjg(zb)/(eps_t(2)-c)
    zvec(1,2) = 1d0/sqrt(1d0+abs(zy)**2)
    zvec(2,2) = zy/sqrt(1d0+abs(zy)**2)

 end if
!write(*,*) zvec
 
end subroutine calc_eig_vec_2x2
!--------------------------------------------------------------------------------------------------

subroutine time_propagation
  use global_variables
  implicit none
integer:: it,ik1,ik2,ik,i
real(8)::Act_t(2)
real(8)::Etot,jt(2),jt_intra2(2)
real(8),allocatable::pop_dist(:,:),berryc(:,:)!pop_dist(nk,2)
!real(8),allocatable :: eps_bk_1d(:,:)
allocate(pop_dist(nk,2))
allocate(berryc(nk,2))
!allocate(eps_bk_1d(2,3*nk))

 call init_laser

 if(myrank==0)then
 open(20,file='total_energy.out')
!   open(21,file='current_intra.out')
   open(22,file='current_total.out')

 end if
 




!do i=0,p
 zpsi=zpsi_gs
do it=0,nt
    !zpsi(t+dt)=exp(-zi*dt*H(t+dt/2)) |zpsi(t)>
    Act_t(:)=0.5d0*(Act(:,it+1)+Act(:,it)) !Act(t+dt/2)=0.5*(Act(t+dt))+Act(t)
    kxt(:)=kx0(:)+Act_t(1)
    kyt(:)=ky0(:)+Act_t(2)

   call dt_evolve
    Act_t(:)=Act(:,it+1)
    kxt(:)=kx0(:)+Act_t(1)
    kyt(:)=ky0(:)+Act_t(2)

!   end do 
call calc_current(jt)
!call calc_intracurrent(jt_intra2)
 call calc_energy(Etot)

  if(myrank==0)then
  write(20,"(2e26.16e3)") dt*it,Etot
!  write(21,"(999e26.16e3)") E0*(i/dble(p)),jt_intra2(:) 
  write(22,"(999e26.16e3)") dt*it,jt(:)

end if
end do 

if(myrank==0)then
    close(20)
 !  close(21)
    close(22)
 end if
   
   call calc_pop_dist(pop_dist)
!    call calc_berry(berryc)

   if(myrank==0)then
  open(21,file='pop_dist_final.out')
  open(23,file='pop_dist_final_a.out')
 !   open(22,file='berry_curvature.out')

   ik=0
do ik1=0,nk1-1
do ik2=0,nk2-1
  ik=ik+1
 write(21,"(999e26.16e3)")kxt(ik),kyt(ik),pop_dist(ik,1:2)
  ! write(22,"(99e26.16e3)")kxt(ik),kyt(ik),berryc(ik,1:2)
end do
write(21,*)
!write(22,*)
end do
close(21)
!close(22)
   ik=0
do ik1=nk1-1, 0, -1
do ik2=nk2-1, 0 ,-1
  ik=ik+1
 write(23,"(999e26.16e3)")kxt(ik),kyt(ik),pop_dist(ik,1:2)
  ! write(22,"(99e26.16e3)")kxt(ik),kyt(ik),berryc(ik,1:2)
end do
write(23,*)
!write(22,*)
end do
close(23)
end if

!call band_1d

! open(20,file='band_1d.out')
!ik=1
!do ik=1,3*nk

 !  write(20,"(999e26.16e3)")dble(ik)/dble(3*nk),eps_bk_1d(1:2,ik)
!end do
!close(20)

end subroutine time_propagation
!-----------------------------------------------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none

  integer:: it,i
  real(8):: tt,xx,theta
real(8):: e_k(2)
  
  allocate(Act(2,-1:nt+1))
  Act=0d0
  e_k(1)=sqrt(3d0)/2d0
  e_k(2)=1d0/2d0
!  do i=0,p
!  theta=i/dble(p)
 theta=0d0*pi/6d0
  do it=0,nt
     tt=dt*it
     xx=tt-0.5d0*Tpulse0
     !yy=tt-0.5d0*Tpulse0
     
     if(abs(xx)<=0.5d0*Tpulse0)then


Act(1,it)=-e_k(1)*(E0)/omega0*(cos(omega0*xx+theta)+1d0/4d0*cos(2d0*omega0*xx+2d0*theta))*cos(pi*xx/Tpulse0)**4
Act(2,it)=-e_k(2)*(E0)/omega0*(cos(omega0*xx+theta)+1d0/4d0*cos(2d0*omega0*xx+2d0*theta))*cos(pi*xx/Tpulse0)**4


  end if
  end do
!  end do

  open(23,file='laser.out')
  do it=0,nt
     tt=dt*it
    write(23,"(999e26.16)")tt,Act(:,it),(Act(:,it+1)-Act(:,it-1))/(2*dt)
  end do
 

 close(23)
  

end subroutine init_laser
!-------------------------------------------------------------------------------------------------

subroutine dt_evolve
  use global_variables
  implicit none
  integer::ik
  complex(8)::zham(2,2),zfk,zvec(2),zhvec(2)
  real(8)::kx_t,ky_t
  complex(8)::zfactor
  integer::iexp

   do ik=nk_s,nk_e
  kx_t=kxt(ik)
  ky_t=kyt(ik)

   zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))

!  ===Taylor expansion ==
  !zpsi(t+dt)=exp(-zi*dt*H(t+dt/2)) |zpsi(t)>
   ! exp(-zi*dt*H(t+dt/2))=1-zi*dtH-0.5d0...
   zfactor=1d0
   zvec(:)=zpsi(:,ik)
   do iexp=1,4
      
      zfactor=zfactor*(-zi*dt)/iexp
      zhvec=matmul(zham,zvec)
      zpsi(:,ik)=zpsi(:,ik)+zfactor*zhvec(:)

      zvec=zhvec

 !  ===Taylor expansion ==
   end do   
end do


end subroutine dt_evolve

!----------------------------------------------------------------------------------------
subroutine calc_energy(Etot)
  use global_variables
  implicit none
  real(8),intent(out):: Etot
 integer::ik
  complex(8)::zham(2,2),zfk,zvec(2),zhvec(2)
  real(8)::kx_t,ky_t
 

  Etot=0d0
  
  do ik=1,nk

     kx_t=kxt(ik)
     ky_t=kyt(ik)

   zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))

   zvec(:)=zpsi(:,ik)
    zhvec(:)=matmul(zham,zvec) ! H|zpsi>

    !Etot=Etot+<zpsi|H|zpsi>
Etot=Etot+sum(conjg(zpsi(:,ik))*zhvec(:))
   
end do
Etot=Etot/nk
end subroutine calc_energy
!------------------------------------------------------------------------------------------------------

subroutine band_1d
  use global_variables
   implicit none
   integer:: ik,ik1,ik2
   complex(8)::zham(2,2),zvec(2,2),zfk
   real(8)::kx_t,ky_t
   real(8):: eps_t(2)
   real(8),allocatable :: eps_bk_1d(:,:)
  ! real(8),intent(out)::eps_bk_1d(2,3*nk)
   allocate(eps_bk_1d(2,3*nk))

!  write(*,*)"ok01"   
   do ik=1,nk
        kx_t=b_vec(1,1)*ik/dble(nk)*2d0/3d0+b_vec(1,2)*ik/dble(nk)*1d0/3d0
        ky_t=b_vec(2,1)*ik/dble(nk)*2d0/3d0+b_vec(2,2)*ik/dble(nk)*1d0/3d0

         zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))
   
   call calc_eig_vec_2x2(zham,zvec,eps_t)
   
   eps_bk_1d(:,ik)=eps_t
     end do

     do ik=1,nk
         kx_t=b_vec(1,1)*(2d0/3d0-ik/dble(nk)*1d0/3d0)+b_vec(1,2)*(1d0/3d0+ik/dble(nk)*1d0/3d0)
         ky_t=b_vec(2,1)*(2d0/3d0-ik/dble(nk)*1d0/3d0)+b_vec(2,2)*(1d0/3d0+ik/dble(nk)*1d0/3d0)
         zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))
   
   call calc_eig_vec_2x2(zham,zvec,eps_t)
  ! zpsi(:,nk+ik)=zvec(:,1)
   eps_bk_1d(:,nk+ik)=eps_t
         
      end do

      do ik=1,nk
         kx_t=b_vec(1,1)*(1-ik/dble(nk))*1d0/3d0+b_vec(1,2)*(1-ik/dble(nk))*2d0/3d0
         ky_t=b_vec(2,1)*(1-ik/dble(nk))*1d0/3d0+b_vec(2,2)*(1-ik/dble(nk))*2d0/3d0
      
   zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))
   
   call calc_eig_vec_2x2(zham,zvec,eps_t)
  ! zpsi(:,2*nk+ik)=zvec(:,1)
   eps_bk_1d(:,2*nk+ik)=eps_t
   
end do


 open(20,file='band_1d.out')

do ik=1,3*nk

   write(20,"(999e26.16e3)")dble(ik)/dble(3*nk),eps_bk_1d(1:2,ik)

! write(20,*)
end do
close(20)

 end subroutine band_1d
!-------------------------------------------------------------------------------------------------------------------------
 subroutine calc_current(jt)
   use global_variables
   implicit none
   real(8),intent(out) :: jt(2)
   real(8)::jt_l(2)
   complex(8)::zJop_x(2,2),zJop_y(2,2),zfk,zvec(2)
   complex(8)::zJvec_x(2),zJvec_y(2)
   real(8):: kx_t,ky_t
   integer::ik
   
   jt_l=0d0
   
   do ik=nk_s,nk_e
   kx_t=kxt(ik)
   ky_t=kyt(ik)

    zJop_x(1,1)=0d0
    zJop_x(2,2)=0d0
    zJop_x(1,2)=-t0_hop*(&
         zi*delta_vec(1,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        +zi*delta_vec(1,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        +zi*delta_vec(1,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))))
    zJop_x(2,1)=conjg(zJop_x(1,2))

    zJop_y(1,1)=0d0
    zJop_y(2,2)=0d0
    zJop_y(1,2)=-t0_hop*(&
        zi*delta_vec(2,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        +zi*delta_vec(2,2)* exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        +zi*delta_vec(2,3)* exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))))
    zJop_y(2,1)=conjg(zJop_y(1,2))

   zvec(:)=zpsi(:,ik)
   zJvec_x=matmul(zJop_x,zvec)
   zJvec_y=matmul(zJop_y,zvec)

   jt_l(1)=jt_l(1)+sum(conjg(zpsi(:,ik))*zJvec_x(:))
   jt_l(2)=jt_l(2)+sum(conjg(zpsi(:,ik))*zJvec_y(:))

end do
call MPI_ALLREDUCE(jt_l,jt,2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
jt=jt/nk
 end subroutine calc_current
!-----------------------------------------------------------------------------------------------
!---------------------------------------------------------------  
subroutine calc_intracurrent(jt_intra2)
  use global_variables
  implicit  none
  real(8),intent(out) :: jt_intra2(2)
  real(8) :: jt_l_intra2(2)
  real(8) :: eps_t(2),square(2),velocity(2,2)
  complex(8) :: zham(2,2), zfk, zphi(2,2)
  complex(8) :: zJop_x(2,2),zJop_y(2,2)
  real(8) :: kx_t, ky_t
  integer :: ik
  real(8)::pop_dist_l(nk,2)
 
 jt_l_intra2 = 0d0

 do ik = nk_s, nk_e

    kx_t = kxt(ik)
    ky_t = kyt(ik)

    zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

    zham(1,1) = eps_b
    zham(1,2) = t0_hop*zfk
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_n  

    call calc_eig_vec_2x2(zham, zphi, eps_t)

    pop_dist_l(ik,1) = abs(sum(conjg(zphi(:,1))*zpsi(:,ik)) )**2
    pop_dist_l(ik,2) = abs(sum(conjg(zphi(:,2))*zpsi(:,ik)) )**2

    zJop_x(1,2) = -t0_hop*( &
      zi*delta_vec(1,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
     +zi*delta_vec(1,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
     +zi*delta_vec(1,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) )

   zJop_y(1,2) = -t0_hop*( &
      zi*delta_vec(2,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
     +zi*delta_vec(2,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
     +zi*delta_vec(2,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) )

   square(1)=-(2d0*eps_t(1)-(eps_b+eps_n))
   square(2)=2d0*eps_t(2)-(eps_b+eps_n)

   velocity(1,1)=((t0_hop**2)/square(1))*&
               (conjg(zfk)*zJop_x(1,2)/(-t0_hop)+conjg(conjg(zfk)*zJop_x(1,2)/(-t0_hop)))
   velocity(2,1)=((t0_hop**2)/square(1))*&
               (conjg(zfk)*zJop_y(1,2)/(-t0_hop)+conjg(conjg(zfk)*zJop_y(1,2)/(-t0_hop)))
   
   velocity(1,2)=-((t0_hop**2)/square(2))*&
               (conjg(zfk)*zJop_x(1,2)/(-t0_hop)+conjg(conjg(zfk)*zJop_x(1,2)/(-t0_hop)))
   velocity(2,2)=-((t0_hop**2)/square(2))*&
               (conjg(zfk)*zJop_y(1,2)/(-t0_hop)+conjg(conjg(zfk)*zJop_y(1,2)/(-t0_hop)))

    jt_l_intra2(1) = jt_l_intra2(1) + pop_dist_l(ik,1)*velocity(1,1) + pop_dist_l(ik,2)*velocity(1,2)
    jt_l_intra2(2) = jt_l_intra2(2) + pop_dist_l(ik,1)*velocity(2,1) + pop_dist_l(ik,2)*velocity(2,2)
  end do

  call MPI_Allreduce(jt_l_intra2, jt_intra2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  jt_intra2 = jt_intra2/nk
end subroutine calc_intracurrent
!---------------------------------------------------------------  

subroutine calc_berry(berryc)
   use global_variables
   implicit none
   real(8),intent(out)::berryc(nk,2)
   real(8),allocatable:: berryc_l(:,:)

   complex(8)::zh_x(2,2),zh_y(2,2)
   complex(8)::zhvec_x1(2),zhvec_y1(2),zhvec_x2(2),zhvec_y2(2)
   real(8):: kx_t,ky_t
   integer::ik
   real(8):: eps_t(2) 
  complex(8)::zham(2,2),zfk,zphi(2,2)
  allocate(berryc_l(nk,2))

   berryc_l=0d0
   
   do ik=nk_s,nk_e
   kx_t=kxt(ik)
   ky_t=kyt(ik)
  
zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))

   call calc_eig_vec_2x2(zham, zphi, eps_t)


    zh_x(1,1)=0d0
    zh_x(2,2)=0d0
    zh_x(1,2)=-t0_hop*(&
         zi*delta_vec(1,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        +zi*delta_vec(1,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        +zi*delta_vec(1,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))))
    zh_x(2,1)=conjg(zh_x(1,2))

    zh_y(1,1)=0d0
    zh_y(2,2)=0d0
    zh_y(1,2)=-t0_hop*(&
        zi*delta_vec(2,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        +zi*delta_vec(2,2)* exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        +zi*delta_vec(2,3)* exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))))
    zh_y(2,1)=conjg(zh_y(1,2))

   
   zhvec_x1=matmul(zh_x,zphi(:,1))
   zhvec_y1=matmul(zh_y,zphi(:,1))
   zhvec_x2=matmul(zh_x,zphi(:,2))
   zhvec_y2=matmul(zh_y,zphi(:,2))

berryc_l(ik,1)=zi*(sum(conjg(zphi(:,1))*zhvec_x2)*sum(conjg(zphi(:,2))*zhvec_y1)&
     -sum(conjg(zphi(:,1))*zhvec_y2)*sum(conjg(zphi(:,2))*zhvec_x1))&
     /((eps_t(1)-eps_t(2))**2)
   
berryc_l(ik,2)=zi*(sum(conjg(zphi(:,2))*zhvec_x1)*sum(conjg(zphi(:,1))*zhvec_y2)&
     -sum(conjg(zphi(:,2))*zhvec_y1)*sum(conjg(zphi(:,1))*zhvec_x2))&
     /((eps_t(1)-eps_t(2))**2)
  
   
end do
call MPI_ALLREDUCE(berryc_l,berryc,2*nk,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

 end subroutine calc_berry


 !-------------------------------------------------------------------------------------------------
 subroutine  calc_pop_dist(pop_dist)
   use global_variables
   implicit none
   real(8),intent(out)::pop_dist(nk,2)
   real(8),allocatable:: pop_dist_l(:,:)

   real(8):: eps_t(2) 
   integer::ik
  complex(8)::zham(2,2),zfk,zvec(2),zphi(2,2)
  real(8)::kx_t,ky_t
  
  allocate(pop_dist_l(nk,2))
  pop_dist_l=0d0
  
  do ik=nk_s,nk_e

     kx_t=kxt(ik)
     ky_t=kyt(ik)

   zfk=exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1)))&
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2)))&
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3)))
   
   zham(1,1)=eps_b
   zham(2,2)=eps_n
   zham(1,2)=t0_hop*zfk
   zham(2,1)=conjg(zham(1,2))

   call calc_eig_vec_2x2(zham, zphi, eps_t)

   pop_dist_l(ik,1) = abs(sum(conjg(zphi(:,1))*zpsi(:,ik)))**2
   pop_dist_l(ik,2) = abs(sum(conjg(zphi(:,2))*zpsi(:,ik)))**2  
   
end do
call MPI_ALLREDUCE(pop_dist_l,pop_dist,2*nk,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
 end subroutine calc_pop_dist
