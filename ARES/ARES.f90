!  ARES.f90 
!
!  FUNCTIONS:
!  ARES - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: ARES
!
!  Корольков Сергей
!
!****************************************************************************

  include "Solver.f90"
	include "STORAGE.f90"
  include "Gas_dyn.f90"
  


program ARES
  USE STORAGE
  USE Government
  implicit none

  call Set_Storage()
  call INIT_Stor()
  call Start(10000)

  call Print_GD()



  pause

end program ARES

