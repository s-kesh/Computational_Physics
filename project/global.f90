module global
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        real, parameter :: a0 = 0.529           ! unit A⁰
        real, parameter :: e = 27.2             ! unit eV hatree energy

        real, parameter :: alpha = 2.
        real, parameter :: gamm = 1.0
        integer, parameter :: ncoord = 3
        integer, parameter :: TotalWalkers = 100
        integer, parameter :: Ntotal = 100000
        integer, parameter :: Nthermal = 10000
        integer, parameter :: minimize_steps = 20
end module
