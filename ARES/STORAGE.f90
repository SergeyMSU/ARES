
module STORAGE
    USE OMP_LIB
  implicit none 

  integer(4), parameter :: n_par = 4                                    ! Число газодинамические параметров параметров
  real(8), parameter :: ggg = (5.0/3.0)
  real(8), parameter :: pi = acos(-1.0_8)        
  real(8), parameter :: alpha = pi/15.0       
  
  TYPE Parameters
    integer :: Nx
    integer :: Ny
    real(8) :: Left
    real(8) :: Right
    real(8) :: Down
    real(8) :: Up
    real(8) :: time_do
    real(8) :: time_posle
    real(8) :: time_all
    real(8) :: dx
    real(8) :: dy
    real(8) :: Volume
    integer (kind=omp_lock_kind) :: time_lock  ! Для openMP

  END TYPE Parameters	

  TYPE Geometry  
  END TYPE Geometry
  
  TYPE GasDyn 
    real(8), allocatable :: GD_Par(:, :, :, :)           ! (4, 2, Nx, Ny) Набор параметров (ro, p, u, v)
  END TYPE GasDyn

  TYPE ATOMS 
    integer N_max        ! Максимальное число атомов в области (какой размер массива)
    integer N_now        ! Максимальное число атомов в области (какой размер массива)
    LOGICAL, allocatable :: AT_be(:)               ! (N_max) Существует ли данный атом в области
    real(8), allocatable :: AT_Par(:, :)           ! (6, N_max) Набор параметров (x, y, Vx, Vy, Vz, I)
  END TYPE ATOMS

  TYPE (Geometry):: GEO1
  TYPE (GasDyn):: GD1
  TYPE (Parameters):: PAR
  TYPE (ATOMS):: AT1

  contains

    subroutine Set_Storage()
        ! Задаём Параметры задачи
        PAR%Left = 0.0
        PAR%Right = 1.0
        PAR%Down = -0.5
        PAR%Up = 0.5
        PAR%Nx = 1500
        PAR%Ny = 1500
        PAR%time_do = 0.00000001
        PAR%time_posle = 0.00000001
        PAR%time_all = 0.0
        PAR%dx = (PAR%Right - PAR%Left)/PAR%Nx
        PAR%dy = (PAR%Up - PAR%Down)/PAR%Ny
        PAR%Volume = PAR%dx * PAR%dy
        call omp_init_lock(PAR%time_lock)

        ! Задаём Газовую динамику
        ALLOCATE(GD1%GD_Par(n_par, 2, PAR%Nx, PAR%Ny))
        GD1%GD_Par = 0.0

        ! Задаём атомы
        AT1%N_max = 1000
        AT1%N_now = 0
        ALLOCATE(AT1%AT_Par(6, AT1%N_max))
        ALLOCATE(AT1%AT_be(AT1%N_max))
        AT1%AT_Par = 0.0
        AT1%AT_be = .False.

    end subroutine Set_Storage

    subroutine Get_koord(n, m, x, y)
        ! Получить координаты центра ячейки по номеру
        INTEGER, INTENT(IN) :: n, m
        real(8), INTENT(OUT) :: x, y

        x = PAR%Left + 1.0 * n/(PAR%Nx + 1) * (PAR%Right - PAR%Left)
        y = PAR%Down + 1.0 * m/(PAR%Ny + 1) * (PAR%Up - PAR%Down)
        return
    end subroutine Get_koord

    subroutine INIT_Stor()

        integer :: i, j
        real(8) :: x, y

        do i = 1, PAR%Nx 
            do j = 1, PAR%Ny
                call Get_koord(i, j, x, y)
                if(y < 0) then
                    GD1%GD_Par(:, 1, i, j) = (/1.0_8, 2.0_8, 0.2_8 * cos(alpha), 0.2_8 * sin(alpha) /)
                    GD1%GD_Par(:, 2, i, j) = GD1%GD_Par(:, 1, i, j)
                else
                    GD1%GD_Par(:, 1, i, j) = (/0.2_8, 2.0_8, 3.0_8 * cos(alpha), 3.0_8 * sin(alpha)/)
                    GD1%GD_Par(:, 2, i, j) = GD1%GD_Par(:, 1, i, j)
                end if
            end do
        end do

    end subroutine INIT_Stor

    subroutine Print_GD()

        integer :: i, j
        real(8) :: x, y, Mach

        open(1, file = 'print_2D.txt')
        write(1,*) "TITLE = 'HP'  VARIABLES = X, Y, RO, P, U, V, Mach"

        do i = 1, PAR%Nx, 3
            do j = 1, PAR%Ny, 3
                call Get_koord(i, j, x, y)
                Mach = sqrt(GD1%GD_Par(3, 1, i, j)**2 + GD1%GD_Par(4, 1, i, j)**2)&
                /sqrt(ggg * GD1%GD_Par(2, 1, i, j)/GD1%GD_Par(1, 1, i, j))
                write(1,*) x, y, GD1%GD_Par(:, 1, i, j), Mach
            end do
        end do

        close(1)


    end subroutine Print_GD
    

end module STORAGE