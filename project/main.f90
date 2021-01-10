program main
        use globaldata
        use twoproton
        implicit none

        real :: S, a
        real, dimension(ncoord * 2, TotalWalkers) :: configuration
        real, dimension(NTotal) :: Elocal, delphibeta, Edelphibeta
        real, dimension(minimize_steps) :: betaarray, earray, evararray
        real, dimension(num_per_node) :: beta_per_node,Energy_per_node,EnergyVar_per_node,S_per_node
        real :: beta, Energy, EnergyVar
        integer :: i, j, loc, seed
        call initmpi

        Nos=num_per_node*numtasks

        call array_allocation


        if(myrank==0)then
                do i=1,Nos
                  S_full_range(i)=0.5 + (i-1)*0.05625
                  !S_full_range(i)=1.0
                end do
        end if
        call mpi_scatter(S_full_range,num_per_node,mpi_real,S_per_node,num_per_node,mpi_real,0,mpi_comm_world,ierror)
        do i = 1, num_per_node
                S = S_per_node(i)
                a = 0.5
                betaarray = 0.4
                call cala(S, a)
                
                do j = 1, minimize_steps
                        call random_seed()
                        call genconfig(configuration)
                        call metropolis(configuration, S, a, betaarray(j), Elocal, delphibeta, Edelphibeta)
                        earray(j) = avg(Elocal(Nthermal: NTotal))
                        evararray(j) = var(Elocal(Nthermal: NTotal))
                        if (j .ne. minimize_steps) &
                                betaarray(j + 1) = min_beta(betaarray(j), earray(j), delphibeta, Edelphibeta)
                end do

                loc = minloc(earray, dim=1)
                beta = betaarray(loc)
                Energy = earray(loc)
                EnergyVar = evararray(loc)

                beta_per_node(i)=beta
                Energy_per_node(i)=Energy
                EnergyVar_per_node(i)=EnergyVar
        end do

        call mpi_gather(beta_per_node,1,mpi_real,beta_gather,1,mpi_real,0,mpi_comm_world,ierror)
        call mpi_gather(Energy_per_node,1,mpi_real,Energy_gather,1,mpi_real,0,mpi_comm_world,ierror)
        call mpi_gather(EnergyVar_per_node,1,mpi_real,EnergyVar_gather,1,mpi_real,0,mpi_comm_world,ierror)

        if(myrank==0)then
                do i=1,Nos
                        write (*, '(A, A, *(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
                          "minimum", TAB, S_full_range(i), TAB, beta_gather(i), TAB, Energy_gather(i), TAB, &
                          EnergyVar_gather(i), N_LINE
                end do
        end if
        call mpi_finalize(ierror)
end program main

subroutine array_allocation
        use globaldata
        implicit none
        if(myrank==0)then
                allocate(S_full_range(Nos))
                allocate(S_gather(Nos))
                allocate(beta_gather(Nos))
                allocate(Energy_gather(Nos))
                allocate(EnergyVar_gather(Nos))
        else
                allocate(S_full_range(0))
                allocate(S_gather(0))
                allocate(beta_gather(0))
                allocate(Energy_gather(0))
                allocate(EnergyVar_gather(0))
        end if
end subroutine

subroutine initmpi
        use globaldata
        implicit none
        call mpi_init(ierror)
        call mpi_comm_rank(mpi_comm_world,myrank,ierror)
        call mpi_comm_size(mpi_comm_world,numtasks,ierror)
        allocate(stats(mpi_status_size))
end subroutine

