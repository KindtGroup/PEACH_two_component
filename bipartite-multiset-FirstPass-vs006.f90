program NVT_from_bulk
! program to generate
!   1. mean concentrations
!   2. equilibrium association constants
!   3. a nucleation free energy profile for a user defined concentration
! Use: bipartite-multiset-inverse.exe < multiset.in > multiset.out
!----------------------------------------
! The anticipated format for multiset.in is:
! #ofDataSets
! #ofIterations alphaValue concentrationForFreeEnergyProfile
! #ofDataPointsFromTraj1 #ofAMoleculesTotal #ofBMoleculesTotal Volume
! #ofAMolecules  #ofBMolecules  #_frequencyAverage_of_NANB_cluster  STDV_in_#frequencyAverage_of_NANB_cluster
! ...
! #ofDataPointsFromTraj1 #ofAMoleculesTotal #ofBMoleculesTotal Volume
! #ofAMolecules  #ofBMolecules  #_frequencyAverage_of_NANB_cluster  STDV_in_#frequencyAverage_of_NANB_cluster
! ...
! etc. for all data sets
!----------------EXAMPLE------------------
! 2                  (NUMBER of DATA SETS)
! 50 0.1 0.05        (FITTING CRITERIA)
! 3 1 1 10.0         (FIRST DATA SET INFORMATION: NLinesData NMolec1 NMolec2 Vol)
! 1 0 1.5  0.003     (FIRST LINE OF DATA: ClusterNMolec1 ClusterNMolec2 ClusterFreq ClusterFreqSTDV)     
! 0 1 1.5  0.003
! 1 1 0.5  0.0006
! 8 2 2 20.0         (SECOND DATA SET INFORMATION)
! 1 0 ...
! 0 1 ...
! 2 0 ...
! 0 2 ...
! 1 1 ...
! 2 1 ...
! 1 2 ...
! 2 2 ...
!----------------------------------------


!----------------------------------------
! DEFINING VARIABLE AND ARRAY TYPES

implicit none ! Cancels default naming conventions
integer, parameter:: Nmax=160
integer, parameter:: Nsetmax=50
       ! to go beyond Nmax=170, will need to 
       ! do something about overflow of factorials

character(8) :: date
character(10) :: time
character(5) :: zone, iterString
character(100) :: startTime,finishTime
real*8 :: startTime_cpu, finishTime_cpu

! Write format, e.g. write(*,FMT1) instead of write (*,*)
character(*), parameter:: FMT1 = "(A1,1x,2(I3,1x),2(1x,E12.5E3))"
character(*), parameter:: FMT2 = "(2(1x,I3),8(1x,E12.5E3))"

real*8:: conc(0:Nmax,0:Nmax,Nsetmax)

real*8:: aveN(0:Nmax,0:Nmax,Nsetmax)
real*8:: aveNout(0:Nmax,0:Nmax,Nsetmax)
real*8:: stdvN(0:Nmax,0:Nmax,Nsetmax)
real*8:: NrelW(0:Nmax,0:Nmax,Nsetmax)
real*8:: weightSum(0:Nmax,0:Nmax)

real*8:: totconc(Nsetmax)
real*8:: V(Nsetmax)
real*8:: Kapp(0:Nmax,0:Nmax,Nsetmax)
real*8:: Kactual(0:Nmax,0:Nmax)
real*8:: KactualCollect(0:Nmax,0:Nmax,0:22+1)
real*8:: alpha
real*8:: nc1

integer:: iset,Nsets,niter,nsamp(0:Nmax,0:Nmax)
integer:: Nsim(Nsetmax,2)
            ! Contains the the total number of each molecule type (Here we 
            ! have limited the number of molecule types to 2 but could be 
            ! expanded should the need arise.)
integer:: nData(Nsetmax) ! Total number of data entries per data set.
integer:: NsimMax(2) ! Has to have two entries for the two molec types
integer:: i,j,k

integer:: FN_LOG, FN_DS

! Keeps track of the monomers in a set that have a conc greater than 0.
logical:: ci_gt_0(0:Nmax,0:Nmax,Nsetmax)
logical:: nsamp_gt_0(0:Nmax,0:Nmax)

100 format (20d8.1)

!----------------------------------------
! GETTING THE DATE AND TIME
call date_and_time(DATE=date,TIME=time,ZONE=zone)
startTime="# Date: "//date(1:4)//"."//date(5:6)//"."//date(7:8)//&
"\n"//" # Time: "//time(1:2)//":"//time(3:4)//":"//time(5:10)//&
" (UTC"//zone(1:3)//":"//zone(4:5)//")"
call cpu_time(startTime_cpu)

write (*,*) startTime
!write(FN_LOG,*) startTime

nsamp_gt_0=.FALSE.
ci_gt_0=.FALSE.
NsimMax=0
totconc=0.0d0
conc=0.0d0
Kapp=0.0d0
Kactual=0.0d0
aveNout=0.0d0
weightSum=0.0d0


!----------------------------------------
! READING DATA

read (*,*) Nsets 
write (*,*) "Number of data sets to be analyzed is ", Nsets
if (Nsets.gt. Nsetmax) then
  write (0,*) "WARNING: Total number of data sets (",Nsets,&
  ") exceeds the maximum (",Nsetmax,")."
  stop
endif
read (*,*) niter, alpha, nc1
! niter: max number of iterations to find best fit equilibrium constants
! alpha: exponent to adjust equilibrium constant at each iteration
! nc1: monomer concentration for DeltaG print in output. This script 
! assumes that the output concentrations of the two components is equal.
! Concentration is given to the program in molec/nm**3.
write (*,*) "Number of iterations is ",niter
write (*,*) "Alpha value for fitting is ", alpha
write (*,*) "Monomer concentration for output is ",nc1," molec/nm^3.\n\n"

do iset=1,nsets
  write (*,*) "iset =", iset
  read (*,*) nData(iset),Nsim(iset,1),Nsim(iset,2),V(iset)
  ! nData: the number of data points in that data set
  ! Nsim: total number of the different molecule types (0 and 1)
  ! V: volume
  
  if (NsimMax(1).lt.Nsim(iset,1)) then 
    NsimMax(1)=Nsim(iset,1)
  endif
  if (NsimMax(2).lt.Nsim(iset,2)) then
    NsimMax(2)=Nsim(iset,2)
  endif
  if (Nsim(iset,1).gt.Nmax) then
    write (0,*) "WARNING: Nsim(",iset,",1) (",Nsim(iset,1),&
    ") is greater than Nmax (",Nmax,").\n    ",&
    "Adjust the Nmax parameter to equal Nsim."
    stop
  endif
  if (Nsim(iset,2).gt.Nmax) then
    write (0,*) "WARNING: Nsim(",iset,",2) (",Nsim(iset,2),&
    ") is greater than Nmax (",Nmax,").\n    ",&
    "Adjust the Nmax parameter to equal Nsim."
    stop
  endif
  
  ! READING MONOMER CONCENTRATIONS FIRST and asserting that it is a
  ! non-zero value. It is ASSUMED that the first two lines of the 
  ! data set pertain to the monomer concentrations. This can be achieved
  ! by sorting the data for each simulation by the total number of 
  ! members in a clusters.
  do k=1,2
    read (*,*) i,j,aveN(i,j,iset),stdvN(i,j,iset)
    if (aveN(i,j,iset).gt.0.0d0) then 
      ci_gt_0(i,j,iset)=.TRUE. ! bool of whether conc > 0 for that i-mer
      nsamp_gt_0(i,j)=.TRUE.
      nsamp(i,j)=nsamp(i,j)+1
    endif
    NrelW(i,j,iset)=aveN(i,j,iset)/stdvN(i,j,iset)
    weightSum(i,j)=weightSum(i,j)+NrelW(i,j,iset)
    conc(i,j,iset)=aveN(i,j,iset)/V(iset)
    totconc(iset)=totconc(iset)+conc(i,j,iset)
    ! apparent equilibrium constant from simulation
    Kapp(i,j,iset)=conc(i,j,iset)/(conc(1,0,iset)**i*conc(0,1,iset)**j)
  enddo
  if (conc(1,0,iset).le.0.0d0) then
    write (0,*) "WARNING: Concentration of free monomer must be > 0.0d0.",&
    "\n    The dataset no. ",iset," has a monomer concentration for molecule",&
    "\n    group 0 of ",conc(1,0,iset)
    stop
  endif
  if (conc(0,1,iset).le.0.0d0) then
    write (0,*) "WARNING: Concentration of free monomer must be > 0.0d0.",&
    "\n    The dataset no. ",iset," has a monomer concentration for molecule",&
    "\n    group 1 of ",conc(0,1,iset)
    stop
  endif

  
  ! READING IN THE REST OF THE I,J-MERS
  do k=3,nData(iset)
    read (*,*) i,j,aveN(i,j,iset),stdvN(i,j,iset)
    if (aveN(i,j,iset).gt.0.0d0) then 
      ci_gt_0(i,j,iset)=.TRUE. ! bool of whether conc > 0 for that i-mer
      nsamp_gt_0(i,j)=.TRUE.
      nsamp(i,j)=nsamp(i,j)+1
    endif
    NrelW(i,j,iset)=aveN(i,j,iset)/stdvN(i,j,iset)
    weightSum(i,j)=weightSum(i,j)+NrelW(i,j,iset)
    ! input bulk system concentration of i-mers
    conc(i,j,iset)=aveN(i,j,iset)/V(iset)
    totconc(iset)=totconc(iset)+(i+j)*conc(i,j,iset)
    ! apparent equilibrium constant from simulation
    Kapp(i,j,iset)=conc(i,j,iset)/(conc(1,0,iset)**i*conc(0,1,iset)**j)
  enddo
  
  write (*,*) "# total monomer concentration is:",totconc(iset)
  write (*,*) "# volume is (",Nsim(iset,1),"+",Nsim(iset,2),")/",totconc(iset)," = ",&
  (Nsim(iset,1)+Nsim(iset,2))/totconc(iset)," =? ",V(iset),"\n"
enddo  ! finish reading data sets

!----------------------------------------
! GENERATE Kactual USING A WEIGHTED GEOMETRIC MEAN OF Kapp
! This is done using a weighted geometric mean of the Kapp values that have been
! read in.
! <x>g = (Sum_i^n x_i^(w_i))^(1/(Sum_i^n w_i))
! where x_i are the i data points and w_i are the respective weightings which in
! this case we have based off of the standard deviation of the individual Kapp
! values which in turn is calculated from the standard devation of the <N_i>.
!----------------------------------------
! 2015.12.14
! In the calculation of Kapp, we must recall that the concentration of the 
! current cluster size is used in tandem with the concentration of the monomers
! such that: Kapp(i) = conc(i)/conc(1)**i
!
! There are two ways that we concidered formulating the weighting parameter
! for each Kapp in the geometric mean:
!       ---------------------------------
!       1) Use the relative weight of the standard deviation in the cluster
!       frequency such that each entry to the geometric mean calculated from
!       Kapp is weighted relative to the data quality of that cluster sizes.
!       The weighting is thus the standard deviation in the cluster frequency
!       divided by said cluster frequency for cluster size i: 
!       NrelW(i)=Nstdv(i)/Navg(i)       (w(i)=1/NrelW(i))
!       REVISION: now NrelW(i)=Navg(i)/Nstdv(i) (w(i)=NrelW(i))
!       ---------------------------------
!       2) Use the relative weight based on the two standard deviations
!       (relative to their respective cluster frequencies) such that for a
!       cluster size i:
!       KrelW(i)=NrelW(i)/(NrelW(1)**i) (w(i)=1/KrelW(i))
!       ---------------------------------
! The idea behind the weighting scheme is to weight the contribution from a 
! data point to the geometric mean for a particular cluster size based on the 
! quality of that data point. 
! In method (2), due to raising the relative
! weighting of the monomer's contribution to the cluster size, this weighting
! tends to over emphasize the monomer's relative standard deviation and as such
! leads to weighting the data sets in their entirity based on the quality of a 
! single data point (the monomer) rather than the data point in question. The
! geometric mean from said weighting shows an unequivocal preference for a
! single data set from i=1-30 and then a sharp drop down to the remaining data
! sets for i>30.
! Method (1) works better in the sense that each data point has an individual
! weight (but does not account for the monomer concentrations data quality in
! the generation of Kapp) and returns a curve that is rather similar to the
! unweighted geometric mean.

do i=0,Nmax
  do j=0,Nmax
      if (Kactual(i,j).eq.0.0d0) then ! that is, if it hasn't been pre-set
        if (weightSum(i,j).ne.0.0d0) then
          Kactual(i,j) = 0.0d0
          do iset=1,Nsets
            if (ci_gt_0(i,j,iset)) then
              Kactual(i,j)=Kactual(i,j)+NrelW(i,j,iset)*log(Kapp(i,j,iset))
            endif
          enddo
          Kactual(i,j)=exp(Kactual(i,j)/weightSum(i,j))
        else
          Kactual(i,j)=0.0d0
        endif
      endif  
  enddo
enddo
KactualCollect(0:Nmax,0:Nmax,0)=Kactual(0:Nmax,0:Nmax)

!----------------------------------------
! STARTING ITERATIONS OF FITTING
call fittingRoutine (Nmax, NsetMax, startTime, &
niter, alpha, Nsets, aveN, aveNout, nsamp_gt_0, &
Nsim, NsimMax, ci_gt_0, NrelW, weightSum, &
V, Kapp, Kactual, KactualCollect, nc1)

!----------------------------------------
! PRINTING OUTPUT
do iset=1,nsets
  write (iterString,"(I0.5)") iset
  open(unit=FN_DS,file='BPMS-DataSet'//iterString//'-OUTPUT.xvg',&
  form='FORMATTED',status='REPLACE',action='WRITE')
  write(FN_DS,*) startTime
  write(FN_DS,*) "# DATASET IDENTIFYING INFORMATION: "
  write(FN_DS,*) "#      Multiset output for data set ", iset
  write(FN_DS,*) "#      Volume = ", V(iset)
  write(FN_DS,*) "#      Number of molecule type 1 = ", Nsim(iset,1)
  write(FN_DS,*) "#      Number of molecule type 2 = ", Nsim(iset,2)
  write(FN_DS,*) "#-------------------------"
  write(FN_DS,*) "# COL | VARIABLE"
  write(FN_DS,*) "#  1  | Number of first molecule type, i"
  write(FN_DS,*) "#  2  | Number of second molecule type, j"
  write(FN_DS,*) "#  3  | Concentration (N_i,j/V) in molec/nm^3 produced by the fit "
  write(FN_DS,*) "#  4  | Concentration (N_i,j/V) from the raw input data "
  write(FN_DS,*) "#  5  | K_actual: The equilibrium association constants from the "
  write(FN_DS,*) "#     | global fitting"
  write(FN_DS,*) "#  6  | K_apparent: The equilibrium association constant calculated"
  write(FN_DS,*) "#     | using the law of mass action with the raw input."
  write(FN_DS,*) "#  7  | Delta G based on K_actual and a monomer concentration of"
  write(FN_DS,*) "#     | ",nc1," molec/nm^3"
  write(FN_DS,*) "#  8  | Data point's weight in the global fitting: "
  write(FN_DS,*) "#     | NrelW(i,j,iset)/weightSum(i,j)."
  write(FN_DS,*) "#-------------------------"
  do i=0,Nsim(iset,1)
    do j=0,Nsim(iset,2)
      if (ci_gt_0(i,j,iset)) then
        write (FN_DS,FMT2) i,j,aveNout(i,j,iset)/V(iset),&
        aveN(i,j,iset)/V(iset),Kactual(i,j),&
        Kapp(i,j,iset),-log(Kactual(i,j))-(i+j-1)*log(nc1),&
        NrelW(i,j,iset)/weightSum(i,j)
      endif
    enddo
  enddo
  close(UNIT=FN_DS)
enddo

!open(UNIT=13,FILE="BPMS-Kactual.xvg",form='FORMATTED',status='REPLACE',action='WRITE')
!write (13,*) startTime
!write (13,*) "# This file contains the equilibrium association constants from "
!write (13,*) "# the final iteration of the global fit."
!write (13,*) "#-------------------------"
!write (13,*) "# COL | VARIABLES"
!write (13,*) "#  1  | Number of first molecule type, i"
!write (13,*) "#  2  | Number of second molecule type, j"
!write (13,*) "#  3  | Association constant after fitting"
!write (13,*) "#  4  | Corresponding Delta_G for a monomer concentration (nc1) of "
!write (13,*) "#     | ",nc1," molec/nm^3: -log(Kactual(i,j))-(i+j-1)*log(nc1)"
!write (13,*) "#-------------------------"
!do i=0,NsimMax(1)
!  do j=0, NsimMax(2)
!    if (nsamp_gt_0(i,j)) then
!      write (13,*) i,j,Kactual(i,j),-log(Kactual(i,j))-(i+j-1)*log(nc1)
!    endif
!  enddo
!enddo
!close(UNIT=13)

!----------------------------------------
! Calculating the time elapsed
call cpu_time(finishTime_cpu)
write (*,*) "Time elapsed = ",(finishTime_cpu-startTime_cpu)/(60.0*60.0)," hrs"

end  ! end main

!----------------------------------------
!--------------SUBROUTINES---------------
!----------------------------------------

subroutine factorial_init (facinv,Nmax)!(factorial, facinv, Nmax)
  implicit none
  !real*8 :: logFacSum(0:Nmax),logFactorial(0:Nmax),logFacInv(0:Nmax),j
  real*8:: factorial(0:Nmax),facinv(0:Nmax)
  integer:: Nmax,i
  ! initialize factorial array
  ! note: floating-point error for Nmax > 170
  ! 2015.12.11 Experimentation using the Exp(Log()) trick has shown that the
  ! sum of log(i)'s remains within the floating-point range of sizes. Once the
  ! exponent of that sum is taken however, we incur the same floating-point 
  ! error for Nmax > 170. This seems non-negotiable unless there's a way to do
  ! all the rest of the calculations using the Exp(log()) format and only take 
  ! the exponent at the very end. 
  factorial(0)=1.0d0
  facinv(0)=1.0d0
  do i=1,Nmax
    factorial(i)=factorial(i-1)*i
    facinv(i)=1.0d0/factorial(i)
  enddo 
  return
end

subroutine bipartite_deriv(n1,n2,ntot,facinv,mat_b,nk)
    ! INPUT:
    ! n1=Nsim(iset,1),n2=Nsim(iset,2)
    ! ntot=Nsim(iset,1)+Nsim(iset,2)
    ! mat_b=qmatrix
    ! facinv(0:ntot)
    ! Note: In upper level, initialize qmatrix(0:n1+1,0:n2+1)=0 
    ! and then give qmatrix input values.
    ! (+1 is to protect code from overflow of data.)
    ! call bipartite_deriv(n1,n2,ntot,facinv(0:ntot),qmatrix(0:n1+1,0:n2+1),nk(0:n1+1,0:n2+1))
    integer,intent(in):: n1,n2,ntot!,Nmax?
    real*8,dimension(0:ntot),intent(in) :: facinv
    real*8,dimension(0:n1+1,0:n2+1),intent(in) :: mat_b
    real*8,dimension(0:n1+1,0:n2+1),intent(out):: nk
    real*8 :: Qn
    integer :: i,j,k1,k2
    real*8,dimension(0:n1+1,0:n2+1) :: mat_a,mat_g,mat_k,mat_c
    real*8,dimension(0:n1+1,0:n2+1) :: mat_x,mat_y,mat_m
    real*8,dimension(0:n1+1,0:n2+1) :: mat_o,C
    real*8,dimension(0:n1+1,0:n2+1) :: deriv_b,deriv_d,deriv_x

    C=0
    Qn=0
    nk=0
    mat_a=0
    mat_c=0
    mat_g=0
    mat_k=0
    mat_x=0
    mat_y=0
    mat_m=0
    mat_o=0
    deriv_x=0
    deriv_b=0
    deriv_d=0
    
    do i=0,n1
        do j=0,n2
            deriv_b(i,j)=mat_b(i+1,j)*(i+1)
        enddo
    enddo
    
    do i=0,n1
        do j=0,n2
            deriv_d(i,j)=mat_b(i,j+1)*(j+1)
        enddo
    enddo
    
    ! Extract Q(N1,N2,V,T) from Xi(lambda1,lambda2,V,T)
    mat_a(0,0)=mat_b(0,0)
    do j=0,n2
        deriv_x(0,j)=mat_b(0,j+1)*(j+1)
    enddo
    mat_a(0,1)=deriv_x(0,0)
    
    mat_x=0
    mat_x=deriv_x
    do j=2,n2
        call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
        mat_a(0,j)=mat_y(0,0)*facinv(j)
        mat_x=0
        mat_x=mat_y
    enddo
    
    mat_m=0
    do j=0,n2
        do i=0,n1
            mat_m(i,j)=mat_b(i+1,j)*(i+1)
        enddo
    enddo
    mat_a(1,0)=mat_m(0,0)
    mat_x=0
    mat_x=mat_m
    do j=1,n2
        call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
        mat_a(1,j)=mat_y(0,0)*facinv(j)
        mat_x=0
        mat_x=mat_y
    enddo

    do i=2,n1
        mat_x=0
        call calc1(n1,n2,deriv_b(0:n1+1,0:n2+1),mat_m(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1))
        mat_a(i,0)=mat_x(0,0)*facinv(i) !! Corrected error here 2017.02.10
        mat_c=mat_x
        do j=1,n2
            mat_y=0
            call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
            mat_a(i,j)=mat_y(0,0)*facinv(i)*facinv(j)
            mat_x=0
            mat_x=mat_y
        enddo
        mat_m=mat_c
        mat_c=0
    enddo
    Qn=mat_y(0,0)*facinv(n1)*facinv(n2)
    
    !Calculate <n_i,j>
    do i=0,n1
        do j=0,n2
            k1=n1-i
            k2=n2-j
            C(k1,k2)=mat_b(k1,k2)*mat_a(i,j)
            nk(k1,k2)=C(k1,k2)/Qn
        enddo
    enddo
end


subroutine calc1(n1,n2,deriv_b,mat_m,mat_x)
    implicit none
    integer,intent(in):: n1,n2
    real*8,dimension(0:n1+1,0:n2+1),intent(in) :: mat_m,deriv_b
    real*8,dimension(0:n1+1,0:n2+1),intent(out)::mat_x
    integer :: r,s,m,j,p,l,k
    real*8,dimension(0:n1+1,0:n2+1) :: deriv_u,mat_u
    real*8,dimension(0:n1*2,0:n2*2) :: mat_v
    real*8,dimension(0:n1+1,0:n2+1) :: mat_w

    mat_x=0
    mat_u=mat_m
    deriv_u=0
    mat_v=0
    mat_w=0
    r=0
    s=0

    do m=0,n2
        do j=0,n1
            deriv_u(j,m)=mat_u(j+1,m)*(j+1)
        enddo
    enddo
    do j=0,n1
        do p=0,n1
            do l=0,n2
                do k=0,n2
                    r=j+p
                    s=l+k
                    mat_v(r,s)=mat_v(r,s)+mat_u(j,l)*deriv_b(p,k)
                enddo
            enddo
        enddo
    enddo
    do l=0,n1
        do m=0,n2
            mat_w(l,m)=deriv_u(l,m)+mat_v(l,m)
        enddo
    enddo
    mat_x=mat_w
end

subroutine calc2(n1,n2,deriv_d,mat_n,mat_y)
    implicit none
    integer,intent(in):: n1,n2
    real*8,dimension(0:n1+1,0:n2+1),intent(in) :: deriv_d,mat_n
    real*8,dimension(0:n1+1,0:n2+1),intent(out)::mat_y
    integer :: s,j,l,k
    real*8,dimension(0:n1+1,0:n2+1) :: deriv_o
    real*8,dimension(0:n1*2,0:n2*2) :: mat_p
    real*8,dimension(0:n1+1,0:n2+1) :: mat_q,mat_o
    mat_y=0
    mat_o=mat_n
    deriv_o=0
    mat_p=0
    mat_q=0
    s=0

    do j=0,n2
        deriv_o(0,j)=mat_o(0,j+1)*(j+1)
    enddo
    do l=0,n2
        do k=0,n2
            s=l+k
            mat_p(0,s)=mat_p(0,s)+mat_o(0,l)*deriv_d(0,k)
        enddo
    enddo
    do l=0,n2
        mat_q(0,l)=deriv_o(0,l)+mat_p(0,l)
    enddo
    mat_y=mat_q
end

subroutine PrintIterN (iter,NsimMax,startTime,nsets,aveNout,Nmax,Nsetmax)
  ! Print iterations data to file.
  implicit none
  character(100) :: startTime
  character(5) :: iterString ! String representation of an integer for file name
  character(14) :: StringAveN
  character(20+(14+1)*50) :: StringAveNs, StringCat
  
  integer :: nsets, NsimMax(2), Nmax, Nsetmax, FN_ITER
  integer :: iter, iset, i, j ! Counters
  
  real*8 :: aveNout(0:Nmax,0:Nmax,Nsetmax)
  
  FN_ITER=20
  write(iterString,"(I0.5)") iter
  if (iter.gt.0) then
    open(UNIT=FN_ITER,FILE="BPMS-Iter"//iterString//"ClusterDistribFit.xvg",&
    form='FORMATTED',status='REPLACE',action='WRITE')
  else
    open(UNIT=FN_ITER,FILE="BPMS-Iter"//iterString//"ClusterDistribRawData.xvg",&
    form='FORMATTED',status='REPLACE',action='WRITE')
  endif
  write(FN_ITER,*) startTime
  write(FN_ITER,*) "# Iteration ",iter
  write(FN_ITER,*) "#-----------------------------------------------"
  write(FN_ITER,*) "# COL | DESCRIPTION"
  write(FN_ITER,*) "#  1  | i, number of molecule type 1"
  write(FN_ITER,*) "#  2  | j, number of molecule type 2"
  write(FN_ITER,*) "#  3  | Data set 1: <n_(i,j)>_(fit) "
  write(FN_ITER,*) "#     | The average frequency that cluster size "
  write(FN_ITER,*) "#     | occurs in data set 1 based on the global "
  write(FN_ITER,*) "#     | fit."
  write(FN_ITER,*) "#  4  | Data set 2: <n_(i,j)>_(fit) "
  write(FN_ITER,*) "# ... "
  write(FN_ITER,*) "# N+2 | Data set N: <n_(i,j)>_(fit) "
  write(FN_ITER,*) "#-----------------------------------------------"
  do i=0,NsimMax(1)
    do j=0,NsimMax(2)
      write (StringAveNs,*) i, j
      do iset=1,nsets
        write(StringAveN,"(ES14.7),1X") aveNout(i,j,iset)
        StringCat=trim(StringAveNs)//StringAveN
        StringAveNs=StringCat
      enddo
      write (FN_ITER,*) trim(StringCat)
    enddo
  enddo
  close(UNIT=FN_ITER)
end

subroutine PrintIterKactuals (KactualCollect,startTime,NsimMax,Nmax,PF,nsamp_gt_0,nc1,niter,Kactual)
  ! Prints out the sequence of saved Kactuals that were collected along the 
  ! fitting process to a single file. Also prints the free energy of associations
  ! to a sister file.
  implicit none
  character(100) :: startTime
  character(14) :: StringK, StringG
  character(1+(14+1)*50) :: StringKs, StringGs, StringL
  
  integer :: NsimMax(2), Nmax, niter
  integer :: i, j, k ! Counters
  integer :: PF, FN_K, FN_G, FN_L
  
  real*8 :: KactualCollect(0:Nmax,0:Nmax,0:22+1), DeltaG, nc1, Kactual(0:Nmax,0:Nmax)
  
  logical :: nsamp_gt_0(0:Nmax,0:Nmax)
  
  FN_K=21
  FN_G=22
  FN_L=23
  open(UNIT=FN_L,FILE="BPMS-LastKactual.xvg",&
  form='FORMATTED',status='REPLACE',action='WRITE')
  
  open(UNIT=FN_K,FILE="BPMS-Kactual.xvg",&
  form='FORMATTED',status='REPLACE',action='WRITE')
  write(FN_K,*) startTime
  write(FN_K,*) "# The equilibrium association constants are printed every ",PF
  write(FN_K,*) "# iterations to this file in collumns."
  write(FN_K,*) "#-----------------------------------------------"
  write(FN_K,*) "# COL | DESCRIPTION"
  write(FN_K,*) "#  1  | i, number of molecule type 1"
  write(FN_K,*) "#  2  | j, number of molecule type 2"
  write(FN_K,*) "#  3  | Iteration 0: Kguess(i,j) from a weighted geometric mean of "
  write(FN_K,*) "#     | the equilibrium association constants approximated "
  write(FN_K,*) "#     | via the law of mass action for that cluster size."
  write(FN_K,*) "#  4  | Iteration 1: Kactual(i,j) from the fit at printing junction 1"
  write(FN_K,*) "#  5  | Iteration ",PF*2,": Kactual(i,j) from the fit at printing junction 2"
  write(FN_K,*) "# ... "
  write(FN_K,*) "# 24  | Iteration ",PF*20,": Kactual(i,j) from the last printing junction 20"
  write(FN_K,*) "# 25  | Iteration ",niter,": Kactual(i,j) from the very last iteration of fitting."
  write(FN_K,*) "#-----------------------------------------------"
  
  open(UNIT=FN_G,FILE="BPMS-DeltaG.xvg",&
  form='FORMATTED',status='REPLACE',action='WRITE')
  write(FN_G,*) startTime
  write(FN_G,*) "# The free energy of association based on the equilibrium "
  write(FN_G,*) "# association constants are printed every ",PF," iterations"
  write(FN_G,*) "# to this file in collumns."
  write(FN_G,*) "# The free energy is calculated assuming a concentration nc1 of "
  write(FN_G,*) "# ",nc1," molec/nm^3 using the equation: "
  write(FN_G,*) "#     -log(Kactual(i,j))-(i+j-1)*log(nc1)"
  write(FN_G,*) "# Kactual(i,j) is the equilibrium association constant for a "
  write(FN_G,*) "# cluster of i molecules of type 1 and j molecules of type 2."
  write(FN_G,*) "#-----------------------------------------------"
  write(FN_G,*) "# COL | DESCRIPTION"
  write(FN_G,*) "#  1  | i, number of molecule type 1"
  write(FN_G,*) "#  2  | j, number of molecule type 2"
  write(FN_G,*) "#  3  | Iteration 0: Free energy of association based on "
  write(FN_G,*) "#     | Kguess(i,j) which is calculated as a weighted geometric "
  write(FN_G,*) "#     | mean of the equilibrium association constants approximated "
  write(FN_G,*) "#     | via the law of mass action for that cluster size."
  write(FN_G,*) "#  4  | Iteration 1: Free energy of association based on "
  write(FN_G,*) "#     | Kactual(i,j) from the fit at printing junction 1."
  write(FN_G,*) "#  5  | Iteration ",PF*1,": Free energy of association based on "
  write(FN_G,*) "#     | Kactual(i,j) from the fit at printing junction 2."
  write(FN_G,*) "# ... "
  write(FN_G,*) "# 24  | Iteration ",PF*20,": Free energy of association based on "
  write(FN_G,*) "#     | Kactual(i,j) from the last printing junction 20."
  write(FN_G,*) "# 25  | Iteration ",niter,": Free energy of association for the "
  write(FN_G,*) "#     | very last iteration of fitting."
  write(FN_G,*) "#-----------------------------------------------"
  
  do i=0,NsimMax(1)
    do j=0,NsimMax(2)
      if (nsamp_gt_0(i,j)) then
        write (StringKs,*) i, j
        write (StringGs,*) i, j
        write (StringL,*) i, j, Kactual(i,j)
        write (FN_L,*) trim(StringL)
        do k=0,22+1
          if (KactualCollect(i,j,k).gt.0) then
            write(StringK,"(ES14.7),1X") KactualCollect(i,j,k)
            StringKs=trim(StringKs)//StringK
            DeltaG = -log(KactualCollect(i,j,k))-(i+j-1)*log(nc1)
            write(StringG,"(F14.7,1X)") DeltaG
            StringGs=trim(StringGs)//StringG
          endif
        enddo
        write (FN_K,*) trim(StringKs)
        write (FN_G,*) trim(StringGs)
      endif
    enddo
    write (FN_G,*) 
  enddo
  close(Unit=FN_K)
  close(Unit=FN_G)
  close(Unit=FN_L)
end

subroutine fittingRoutine (Nmax, NsetMax, startTime, &
niter, alpha, Nsets, aveN, aveNout, nsamp_gt_0, &
Nsim, NsimMax, ci_gt_0, NrelW, weightSum, &
V, Kapp, Kactual, KactualCollect, nc1)
! 2015.12.28 : This is the subscript to which the entire fitting scheme was
! moved to.
  implicit none
  
  character(100) :: startTime
  
  integer :: FN_LOG, FN_CONV !! Unit numbers for writing files
  integer :: PF ! Printing frequency
  integer :: Nmax ! Maximum i-mer size
  integer :: NsetMax ! Maximum number of data sets
  integer :: niter ! Number of iterations
  integer :: Nsets ! Number of data sets
  integer :: iter, iset, i, j ! Counters
  integer :: nsamp(0:Nmax,0:Nmax) ! Number of samples for each i-mer
  integer :: Nsim(Nsetmax,2), NsimMax(2)
  integer :: countConverCrit

  real*8 :: aveN(0:Nmax,0:Nmax,Nsetmax), aveNout(0:Nmax,0:Nmax,Nsetmax)
  real*8 :: NrelW(0:Nmax,0:Nmax,Nsetmax), weightSum(0:Nmax,0:Nmax)
  real*8 :: V(Nsetmax)
  real*8 :: Kapp(0:Nmax,0:Nmax,Nsetmax), Kactual(0:Nmax,0:Nmax)
  real*8 :: KactualCollect(0:Nmax,0:Nmax,0:22+1)
  
  real*8 :: alpha, Kapp_current, nc1
  real*8 :: facinv(0:Nmax), q(0:Nmax,0:Nmax)
  real*8 :: Kfac, aveNoutC(0:Nmax,0:Nmax)
  real*8 :: converVar
  real*8 :: converCritGM(0:Nmax,0:Nmax),converCritTotGM

  integer :: NsimC(2),printCount
  integer :: n1,n2,ntot

  ! Bools that control whether or not to print in that iteration.
  ! Evaluated only once at the begining of the iteration.
  logical :: printing
  logical :: ci_gt_0(0:Nmax,0:Nmax,Nsetmax) 
  logical :: nsamp_gt_0(0:Nmax,0:Nmax)
 
  FN_CONV=12
  aveNoutC=0.0d0
  printCount=0

  ! Printing frequency
  if (niter.gt.20) then
    PF=niter/20
  else
    PF=1
  endif
  
  !----------------------------------------
  ! OPENING FILES AS NEEDED
  ! Note that the standard unit numbers for fortran are:
  ! Standard Error = 0 : Used to print error messages to the screen.
  ! Standard In = 5    : Used to read in data from the keyboard.
  ! Standard Out = 6   : Used to print general output to the screen.
  open(unit=FN_LOG,file='BPMS-LOG.log',form='FORMATTED',status='REPLACE',action='WRITE')
  open(UNIT=FN_CONV,FILE="BPMS-ConvergenceCriteria.xvg",&
  form='FORMATTED', status='REPLACE',action='WRITE')
  write(FN_CONV,*) startTime
  write(FN_CONV,*) "# This file contains a plottable progression for the "
  write(FN_CONV,*) "# convergence criteria, primarily the result of sum: "
  write(FN_CONV,*) "#   ((n_model(i,iset)-n(i,iset))**2/n_stdv(i,iset)**2)"
  write(FN_CONV,*) "# over all the data. The desired trend should be to see"
  write(FN_CONV,*) "# the value of this go down and the iteration step goes"
  write(FN_CONV,*) "# up."
  write(FN_CONV,*) "# COL | VARIABLE "
  write(FN_CONV,*) "#  1  | iteration step"
  write(FN_CONV,*) "#  2  | sum_(i) ( prod_(iset) ("
  write(FN_CONV,*) "#     | ((n_model(i,iset)-n(i,iset))**2/n_stdv"
  write(FN_CONV,*) "#     | (i,iset)**2)**(NrelW(i,iset)/weightSum(i))"
  write(FN_CONV,*) "#     | ) )"
  write(FN_CONV,*) "# "
  
  call PrintIterN (0,NsimMax,startTime,nsets,aveN,Nmax,Nsetmax)
  
  ! Defining the factorials and inverse factorials using the module/function
  ! factorial_init
  call factorial_init (facinv,Nmax)  
  
  do iter=1,niter
    ! To print or not to print
    printing=((mod(iter,PF).eq.0).or.(iter.eq.1))
    if (printing) then
      write (*,*) "Iteration ",iter
      write(FN_LOG,*) "#---------------------"
      write(FN_LOG,*) "# Iteration = ", iter
      write(FN_LOG,*) "#---------------------"
      write(FN_LOG,*) "# i, j, converCritGM(i,j)"
    endif

    do iset=1,Nsets
      NsimC=Nsim(iset,:) ! total number of monomers in this run
      q(0,0)=1.0d0
      q(0,1)=1.0d0
      q(1,0)=1.0d0
      do i=0,Nmax
        do j=0,Nmax
        ! generate "q" (actually scaled to V) array from current best
        ! guess Kactual
        !----------------------------------------
        ! INITIAL GUESS OF q(i) FROM APPARENT K - SINGLE COMPONENT
        ! First guess at finding single-aggregate partition functions relative
        ! to monomer. (For a single data set, this will set Kactual = Kapp)
        !----------------------------------------
        ! Note: q0(i) is the standard state partition function for an i-mer.
        !       V0 is the standard state volume, assumed to be equal to 1.
        !       K(i) is the equilibrium constant for the i-mer and equals
        !               q0(i)/q0(1)**i
        ! q(i)=q0(i)*(V/V0)
        !     =K(i)*(q0(1)**i/q0(i))*q0(i)*(V/V0)  
        !               !multiplying by 1=K(i)*K(i)**-1
        !     =K(i)*q0(1)**i*(V/V0)
        !     =K(i)*q(1)**i*(V0/V)**i*(V/V0)       !q0(1)=q(1)*(V0/V)
        !     =K(i)*q(1)**i*(V0/V)**(i-1)
        !     =K(i)*q(1)**i*V**(1-i)               !V0=1
        !     =K(i)*V**(1-i)                       !q(1)=1.0d0
        !----------------------------------------
        ! INITIAL GUESS OF q(a,b) FROM APPARENT K - TWO COMPONENT
        !----------------------------------------
        ! Note: q0(a,b) is the standard state partition function for an a,b-mer.
        !       V0 is the standard state volume, assumed to be equal to 1 nm**3.
        !       K(a,b) is the equilibrium constant for the a,b-mer and equals
        !               q0(a,b)/(q0(1,0)**a*q0(0,1)**b)
        ! q(a,b)=q0(a,b)*(V/V0)
        !       =K(a,b)*(q0(1,0)**a*q0(0,1)**b)/q0(a,b)*q0(a,b)*(V/V0)
        !               !multiplying by 1=K(a,b)*K(a,b)**-1
        !       =K(a,b)*(q0(1,0)**a*q0(0,1)**b)*(V/V0)
        !               !q0=q*(V0/V) substituted in for q0(1,0) and q0(0,1)
        !       =K(a,b)*(q(1,0)*(V0/V))**a*(q(0,1)*(V0/V))**b**(V/V0)
        !       =K(a,b)*q(1,0)**a*q(0,1)**b*(V0/V)**(a+b-1)
        !               !Assuming the   q(1,0)=1.0d0
        !               !               q(0,1)=1.0d0
        !               !               V0=1 nm**3
        !       =K(a,b)*V**(1-(a+b))
        !----------------------------------------
          if (((i+j).gt.1).and.(ci_gt_0(i,j,iset))) then 
            q(i,j)=Kactual(i,j)*V(iset)**(1-i-j)
          else if ((i+j).gt.1) then
            q(i,j)=0.0d0
          endif 
          
        enddo
      enddo
      n1=Nsim(iset,1)
      n2=Nsim(iset,2)
      ntot=n1+n2
      ! Subroutine forward returns aveNoutC
      
      call bipartite_deriv(n1,n2,ntot,facinv(0:ntot),q(0:n1+1,0:n2+1),&
      aveNoutC(0:n1+1,0:n2+1))
      aveNout(0:n1+1,0:n2+1,iset)=aveNoutC
    enddo
    ! Tracking the convergence of Kactual
    countConverCrit=0
    converCritGM=0.0d0
    converCritTotGM=0.0d0
    do i=0,Nmax
      do j=0,Nmax
        if ((i+j).gt.1) then
          do iset=1,nsets
            if (ci_gt_0(i,j,iset)) then
              ! The following calculation keeps track (in a loose sense) of the
              ! quality of the fitting for each i-mer individually and also the
              ! data set as a whole. Contributions here are again weighted by
              ! NrelW/weightSum.
              ! GEOMETRIC AVERAGE CONVERSION CRITERIA
              converVar=(NrelW(i,j,iset)/weightSum(i,j))*&
              log((aveNout(i,j,iset)-aveN(i,j,iset))**2/aveN(i,j,iset)**2)
              converCritGM(i,j)=converCritGM(i,j)+converVar
            endif
          enddo
          if (converCritGM(i,j).ne.0.0d0) then
            countConverCrit=countConverCrit+1
            converCritGM(i,j)=exp(converCritGM(i,j))
            if (printing) then
              write (FN_LOG,*) i, j, converCritGM(i,j)
            endif
          endif
          converCritTotGM=ConverCritTotGM+converCritGM(i,j)
        endif
      enddo
    enddo
    converCritTotGM=converCritTotGM/dfloat(countConverCrit)
    write(FN_CONV,*) iter, converCritTotGM
    !! END OF CONVERSION CRITERIA CALCULATION
    
    do i=0,Nmax
      do j=0,Nmax
        if ((i+j).gt.1) then
          Kfac=0.0d0
          if (nsamp_gt_0(i,j)) then
            do iset=1,nsets
              if (ci_gt_0(i,j,iset)) then
                ! calculate the apparent association constant from the latest
                ! iteration. Make a new guess based on the ratio of the current 
                ! constant to the true value; each set was originally weighted 
                ! equally:
                !     Kfac=Kfac*(Kapp(i,iset)/Kapp_current)**(alpha/dfloat(nsamp))
                ! but that has been changed using NrelW/weightSum which
                ! should appropriate a greater weight to data points that have a
                ! smaller relative standard deviation.
                Kapp_current=(V(iset)**(i+j-1))*aveNout(i,j,iset)/(aveNout(1,0,iset)**i*aveNout(0,1,iset)**j)
                Kfac=Kfac+log(Kapp(i,j,iset)/Kapp_current)*(alpha*&
                NrelW(i,j,iset)/weightSum(i,j))
              endif
            enddo
            ! Scaling the Kactual value by its corresponding exp(Kfac)
            Kactual(i,j)=Kactual(i,j)*exp(Kfac) 
          endif
        endif  
      enddo
    enddo
    
    if (printing) then
      printCount=printCount+1
      KactualCollect(0:Nmax,0:Nmax,printCount)=Kactual(0:Nmax,0:Nmax)
      call PrintIterN (iter,NsimMax,startTime,nsets,aveNout,Nmax,Nsetmax)
    endif
  enddo
  close(UNIT=FN_CONV)
  close(UNIT=FN_LOG)
  printCount=printCount+1
  KactualCollect(0:Nmax,0:Nmax,printCount)=Kactual(0:Nmax,0:Nmax)
  call PrintIterKactuals (KactualCollect,startTime,NsimMax,Nmax,PF,nsamp_gt_0,nc1,niter, Kactual)
  return
end
