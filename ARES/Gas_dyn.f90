module Government
    USE STORAGE
    USE Solver
    USE OMP_LIB
    implicit none 
    contains

    subroutine Start(step_all)
        integer, INTENT(IN) :: step_all
        integer :: step, now, now2, i, n, m

        now = 2
        now2 = 1

        do step = 1, step_all 
            if(mod(step, 1000) == 0) then
                print*, "step = ", step, "d_time = ", PAR%time_posle, now, now2
            end if
            now = mod(now, 2) + 1      ! Какие переменные сейчас используем
            now2 = mod(now2, 2) + 1    ! Какие переменные вычисляем
            PAR%time_do = PAR%time_posle
            PAR%time_posle = 1000.0
            ! Пробегаемся по всем ячейкам и вычисляем газодинамические параметры на следующем шаге
            !$omp parallel
            !$omp do private(n, m)
            do i = 1, PAR%Nx * PAR%Ny
                n = mod(i - 1, PAR%Nx) + 1   ! Нашли номер ячейки в формате (n, m)
                m = (i - n) / PAR%Nx + 1
                call Culc_Cell(n, m, now, now2)  ! Находим газодинамические параметры в ячейке на следующем шаге
            end do
            !$omp end do
            !$omp end parallel

        end do

    end subroutine Start


    subroutine Culc_Cell(n, m, now, now2)
        ! Полностью вычислияем газодинамические параметры в ячейке на следующем шаге
        ! Внутри вычислияется минимальный шаг по времени
        integer, INTENT(IN) :: n, m, now, now2
        real(8) :: s1(n_par), s2(n_par), s3(n_par), s4(n_par), s5(n_par)
        real(8) :: x, y, dsl, dsp, dsc, qqq1(5), qqq2(5), qqq(5), POTOK(5)
        real(8) :: ro, u, v, p, loctime

        call Get_koord(n, m, x, y)
        s1 = GD1%GD_Par(:, now, n, m)
        POTOK = 0.0
        loctime = 100000.0

        if(n == 1) then
            if(y < 0) then
                s4 = (/1.0_8, 2.0_8, 0.2_8 * cos(alpha), 0.2_8 * sin(alpha) /)         ! Граничные условия
            else
                s4 = (/0.2_8, 2.0_8, 3.0_8 * cos(alpha), 3.0_8 * sin(alpha) /)         ! Граничные условия
            end if
        else
            s4 = GD1%GD_Par(:, now, n - 1, m)
        end if

        if(n == PAR%Nx) then
            s2 = s1         ! мягкие раничные условия
        else
            s2 = GD1%GD_Par(:, now, n + 1, m)
        end if

        if(m == 1) then
            !s3 = s1         ! мягкие раничные условия
            s3 = (/1.0_8, 2.0_8, 0.2_8 * cos(alpha), 0.2_8 * sin(alpha) /)         ! Граничные условия
        else
            s3 = GD1%GD_Par(:, now, n, m - 1)
        end if

        if(m == PAR%Ny) then
            !s5 = s1         ! мягкие раничные условия
            s5 = (/0.2_8, 2.0_8, 3.0_8 * cos(alpha), 3.0_8 * sin(alpha) /)         ! Граничные условия
        else
            s5 = GD1%GD_Par(:, now, n, m + 1)
        end if

        qqq1(1) = s1(1); qqq1(2) = s1(3); qqq1(3) = s1(4); qqq1(4) = 0.0; qqq1(5) = s1(2)

        qqq2(1) = s2(1); qqq2(2) = s2(3); qqq2(3) = s2(4); qqq2(4) = 0.0; qqq2(5) = s2(2)
        call chlld_gd(2, 1.0_8, 0.0_8, 0.0_8, 0.0_8, qqq1, qqq2, dsl, dsp, dsc, qqq)
        POTOK = POTOK + qqq * PAR%dy
        loctime = min(loctime, 0.5 * PAR%dx/max(dabs(dsl), dabs(dsp)))

        qqq2(1) = s3(1); qqq2(2) = s3(3); qqq2(3) = s3(4); qqq2(4) = 0.0; qqq2(5) = s3(2)
        call chlld_gd(2, 0.0_8, -1.0_8, 0.0_8, 0.0_8, qqq1, qqq2, dsl, dsp, dsc, qqq)
        POTOK = POTOK + qqq * PAR%dx
        loctime = min(loctime, 0.5 * PAR%dy/max(dabs(dsl), dabs(dsp)))

        qqq2(1) = s4(1); qqq2(2) = s4(3); qqq2(3) = s4(4); qqq2(4) = 0.0; qqq2(5) = s4(2)
        call chlld_gd(2, -1.0_8, 0.0_8, 0.0_8, 0.0_8, qqq1, qqq2, dsl, dsp, dsc, qqq)
        POTOK = POTOK + qqq * PAR%dy
        loctime = min(loctime, 0.5 * PAR%dx/max(dabs(dsl), dabs(dsp)))

        qqq2(1) = s5(1); qqq2(2) = s5(3); qqq2(3) = s5(4); qqq2(4) = 0.0; qqq2(5) = s5(2)
        call chlld_gd(2, 0.0_8, 1.0_8, 0.0_8, 0.0_8, qqq1, qqq2, dsl, dsp, dsc, qqq)
        POTOK = POTOK + qqq * PAR%dx
        loctime = min(loctime, 0.5 * PAR%dy/max(dabs(dsl), dabs(dsp)))

        if(loctime < PAR%time_posle) then
            call omp_set_lock(PAR%time_lock)
                if(loctime < PAR%time_posle) then
                    PAR%time_posle = loctime
                end if
            call omp_unset_lock(PAR%time_lock)
        end if

        ro = s1(1) - PAR%time_do/PAR%Volume * POTOK(1)
        if(ro <= 0.0) then
            print*, "Ro < 0", s1(1), ro
            print*, x, y, n, m
            STOP
        end if
        u = (s1(1) * s1(3) - PAR%time_do/PAR%Volume * POTOK(2))/ro
        v = (s1(1) * s1(4) - PAR%time_do/PAR%Volume * POTOK(3))/ro
        p = (((s1(2) / (ggg - 1.0) + s1(1) * (s1(3)**2 + s1(4)**2) * 0.5) - (PAR%time_do/PAR%Volume) * POTOK(5)) - &
        0.5 * ro * (u**2 + v**2)) * (ggg - 1.0);
        if(p <= 0.0000001) then
            p = 0.0000001
        end if

        GD1%GD_Par(:, now2, n, m) = (/ro, p, u, v/)

    end subroutine Culc_Cell






end module Government