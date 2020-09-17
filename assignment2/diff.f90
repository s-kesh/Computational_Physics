program diffeqnsolver
        use ode1d, only : odesolve
        use function
        implicit none
        real :: x, y, x_final, step_size

        x = 0
        y = 1
        step_size = 0.001
        x_final = 4

        call odesolve(f, "SE", x, y, x_final, step_size)
        call odesolve(f, "ME", x, y, x_final, step_size)
        call odesolve(f, "RK", x, y, x_final, step_size)

end program diffeqnsolver

