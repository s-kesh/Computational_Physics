module function
        implicit none
        Public f
        contains
                real function f(x,y)
                        implicit none
                        real, intent(in) :: x, y
                        f = -2*y
                end function
end module function
