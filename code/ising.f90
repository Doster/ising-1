! Name:      Jared Doster and Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Monte Carlo Ising Model: Wolff Algorithm and Metropolis Test
! 
!Program Summary (Ising)
!
!  Input:    Hardcode three parameters to the subroutine (temperature of electrons
!            and the length and width of the eletron array).  
!
!  Process:  Creates an array representing the arry of electrons.
!            Calls the subroutine in a nested loop from temp=0 to tempfinal=arbitrary.
!            Each iteration of the outer "temperature" loop returns the magnetizaton for a
!            particular temperature.
!            The iterations of the inner loop will continue until the magnetization returned
!            from the subroutine has converged.
!
!  Output:   Prints magnetization vs. temperature and magnetization vs. time into two data files
!
!
!Subroutine Summary (Mainloop)
!
! Input:     Three parameters from the program
!
! Process:   For a particular temperature, runs a loop involving random selection of electrions
!            and performing a Metropolis Test. Each time that the mainloop is called, the spins
!            of the array are randomly chosen and flipped.
!            Many iterations of this loop cause the calculated magnetization to converge to a range
!            of values.
!
! Output:    Magnetization for a particular temperature for a particular iteration.
!            Output will vary each time that the subroutine is called



program ising

  use model
  use plot
 
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INPUT: Row and column size
  integer,parameter :: SIZE = 32

!! INPUT: Final temperature (Kelvin x 100), final time, stepsizeloops
  integer,parameter :: TEMPFINAL = 350, TIMEFINAL = 100000

!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:SIZE-1,0:SIZE-1), temp, time
  real(8) :: weight(-2:2), spintotal, spinmagavg2
  integer :: wolffspin(0:SIZE-1, 0:SIZE-1)
  real(8) :: wolfftotal, wolffmagavg2
  integer :: swenwangspin(0:SIZE-1, 0:SIZE-1)
  real(8) :: swenwangtotal, swenmagavg2
!  real(8) :: swenenergy, metropenergy, wolffenergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call plot_init()

!!! Open files for writing data !!!                                                                                                                                      
  call opentextfiles

!! Run main loop subroutines (both metropolis and Wolff run simultaneously)
!! In the following nested loops, the outer loop is the temperature loop 
!!  (for plot of magnetization vs. temperature) 
!! and the inner loop is the converging time loop 
!!  (for plot of magnetization vs. time)


  do temp = 150,TEMPFINAL


! Re-initialize Wolff lattice. All values in the lattice must be 1
    wolffspin(:,:) = 1
    wolfftotal = 0
    wolffmagavg2 = 0
    do time = 1, TIMEFINAL, 100
      call wolff(wolffspin, SIZE, temp/100d0)
      wolffmagavg2 = wolffmagavg2 + wolfftotal**2
      wolfftotal = wolfftotal + abs(sum(wolffspin)/(SIZE**2*1d0))
    end do
    wolfftotal = 100d0 * wolfftotal / TIMEFINAL
    wolffmagavg2 = 100d0* wolffmagavg2 / TIMEFINAL

!    wolffenergy = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!    call energy(wolffspin, SIZE, wolffenergy) !!!!!!!!!!!!!!!
    

    
! Re-initialize Swendswon-Wang lattice. All values in the lattice must be 1
    swenwangspin(:,:) = 1
    swenwangtotal = 0
    swenmagavg2 = 0
    do time = 1, TIMEFINAL, 100
      call swenwang(swenwangspin, SIZE, temp/100d0)
      swenmagavg2 = swenmagavg2 + swenwangtotal**2
      swenwangtotal = swenwangtotal + abs(sum(swenwangspin)/(SIZE**2*1d0))
    end do
    swenwangtotal = 100d0 * swenwangtotal / TIMEFINAL
    swenmagavg2 = 100d0* swenmagavg2 / TIMEFINAL

!    swenenergy = 0
!    call energy(swenwangspin, SIZE, swenenergy)



! Re-initialize the Metropolis lattice. All values in the lattice must be 1
    spin(:,:) = 1
    spintotal = 0       
    spinmagavg2 = 0
! Only 5 options for the exponent for metropolis so calculate them once
    weight = [exp(-800d0/temp),exp(-400d0/temp),1d0,exp(400d0/temp),exp(800d0/temp)]
    do time = 0,TIMEFINAL
      call metropolis(spin, SIZE, weight)
      spinmagavg2 = spinmagavg2 + spintotal**2
      spintotal = spintotal + abs(sum(spin)/(SIZE**2*1d0))
! We want time to print only once (choose an arbitrary temperature)
      if (temp == 250) then
        WRITE(15,*) sum(spin)/(SIZE**2*1d0), time
      end if
    end do
    spintotal = spintotal / TIMEFINAL
    spinmagavg2 = spinmagavg2 / TIMEFINAL

!    metropenergy = 0 
!    call energy(spin, SIZE, metropenergy)

     call plot_spin(wolffspin, SIZE, temp/100d0)

! Print magnetization
    WRITE(16,*) spintotal, temp/100d0
    WRITE(17,*) swenwangtotal, temp/100d0
    WRITE(18,*) wolfftotal, temp/100d0

! Calculate magnetic susceptibility and print
    WRITE(19,*) (100d0/temp)*( spinmagavg2 - spintotal**2 ), temp/100d0
    WRITE(20,*) (100d0/temp)*( wolffmagavg2 - wolfftotal**2 ), temp/100d0
    WRITE(21,*) (100d0/temp)*( swenmagavg2 - swenwangtotal**2 ), temp/100d0

! Print energy
!    WRITE(22,*) metropenergy, temp/100d0 
!    WRITE(23,*) wolffenergy, temp/100d0 
!    WRITE(24,*) swenenergy, temp/100d0 

  end do


!!! Close text files !!!                                                                                                                                                 
  call closetextfiles
  call plot_close()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="metrop_mag_time.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_mag_time file not opened properly------------"
    endif

    OPEN(UNIT=16,FILE="metrop_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_mag_temp file not opened properly------------"
    endif
    OPEN(UNIT=17,FILE="swenwang_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, swenwang_mag_temp file not opened properly------------"
    endif
    OPEN(UNIT=18,FILE="wolff_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, wolff_mag_temp file not opened properly------------"
    endif

    OPEN(UNIT=19,FILE="metrop_magsus_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_magsus_temp file not opened properly------------"
    endif
    OPEN(UNIT=20,FILE="swenwang_magsus_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, swenwang_magsus_temp file not opened properly------------"
    endif
    OPEN(UNIT=21,FILE="wolff_magsus_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, wolff_magsus_temp file not opened properly------------"
    endif

!    OPEN(UNIT=22,FILE="metrop_energy.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, metrop_energy file not opened properly------------"
!    endif
!    OPEN(UNIT=23,FILE="wolff_energy.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, wolff_energy file not opened properly------------"
!    endif
!    OPEN(UNIT=24,FILE="swen_energy.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, swen_energy file not opened properly------------"
!    endif

!    OPEN(UNIT=25,FILE="metrop_spec.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, metrop_spec file not opened properly------------"
!    endif
!    OPEN(UNIT=26,FILE="wolff_spec.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, wolff_spec file not opened properly------------"
!    endif
!    OPEN(UNIT=27,FILE="swen_spec.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
!    if (OPEN_STATUS /= 0) then
!       STOP "------------Error, swen_spec file not opened properly------------"
!    endif
end subroutine


subroutine closetextfiles
    CLOSE(UNIT=15)

    CLOSE(UNIT=16)
    CLOSE(UNIT=17)
    CLOSE(UNIT=18)

    CLOSE(UNIT=19)
    CLOSE(UNIT=20)
    CLOSE(UNIT=21)

!    CLOSE(UNIT=22) 
!    CLOSE(UNIT=23) 
!    CLOSE(UNIT=24) 

!    CLOSE(UNIT=25) 
!    CLOSE(UNIT=26) 
!    CLOSE(UNIT=27) 
end subroutine

end program ising
