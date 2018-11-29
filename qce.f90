!=========================================================================================
! Peacemaker -- A Quantum Cluster Equilibrium Code.
!
! Copyright 2004-2006 Barbara Kirchner, University of Bonn
! Copyright 2007-2012 Barbara Kirchner, University of Leipzig
! Copyright 2013-2018 Barbara Kirchner, University of Bonn
!
! This file is part of Peacemaker.
!
! Peacemaker is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Peacemaker is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
!=========================================================================================
! This module implements the QCE algorithm, the core of Peacemaker.
module qce
    use omp_lib
    use kinds
    use cluster
    use constants
    use auxiliary, only: range_t
    use shared_data
    use partition_functions
    implicit none
    private
    !=====================================================================================
    ! Public entities
    public :: qce_prepare
    public :: qce_start
    public :: qce_finalize
    !=====================================================================================
    ! Data type storing reference_data.
    type :: reference_t
        logical :: compare, compare_isobar, compare_density, compare_phase_transition
        real(dp) :: density_weight, isobar_weight, phase_transition_weight
        real(dp) :: phase_transition, density, density_temperature
        real(dp), dimension(:), allocatable :: isobar_temperature, isobar_volume
    end type reference_t
    !=====================================================================================
    ! Instantiation of the data type above.
    type(reference_t):: reference
    !=====================================================================================
    ! Data type representing an isobar. Instantions of this data type are distributed
    ! over threads. Together with global_data, they fully describe a QCE calculation.
    type :: isobar_t
        real(dp) :: amf, bxv, error
        integer, dimension(:), allocatable :: solution
        logical, dimension(:), allocatable :: converged
        real(dp), dimension(:), allocatable :: vol, temp, gibbs
        real(dp), dimension(:, :), allocatable :: populations
        type(pf_t), dimension(:, :), allocatable :: lnq
    end type isobar_t
    !=====================================================================================
    ! Instantiation of the data type specified above. This is the data that qce_main()
    ! works with.
    type(isobar_t):: ib
    type(isobar_t):: best_ib
    !=====================================================================================
    contains
        !=================================================================================
        ! Prepares the isobar array. Use of input data is only legit inside this
        ! subroutine. Unit conversion may be performed here.
        subroutine qce_prepare()
            ! Use of input data is only legit inside this subroutine.
            use input
    
            ! Assign global data.
            global_data%press = pmk_input%pressure*1.0e5_dp ! Conversion to Pa
            global_data%max_deviation = pmk_input%max_deviation
            global_data%vdamp = pmk_input%volume_damping_factor
            global_data%qce_iterations = pmk_input%qce_iterations
            global_data%newton_iterations = pmk_input%newton_iterations
            global_data%amf = pmk_input%amf
            global_data%bxv = pmk_input%bxv
            global_data%temp = pmk_input%temperature
            global_data%progress_bar = pmk_input%progress_bar
    
            allocate(global_data%amf_pure(size(monomer)))
            global_data%amf_pure = pmk_input%amf_pure/avogadro**2
            allocate(global_data%bxv_pure(size(monomer)))
            global_data%bxv_pure = pmk_input%bxv_pure
            allocate(global_data%monomer_amounts(size(monomer)))
            global_data%monomer_amounts = pmk_input%monomer_amounts
            call initialize_conserved_quantities()

            allocate(global_data%degree(size(monomer)))
            call initialize_degree()
    
            global_data%nconverged = 0
    
            ! Assign reference data.
            reference%compare = pmk_input%compare
    
            reference%compare_isobar = pmk_input%compare_isobar
            if (reference%compare_isobar) then
                reference%isobar_weight = pmk_input%ref_isobar_weight
                allocate(reference%isobar_temperature( &
                    size(pmk_input%ref_isobar_temperature)))
                reference%isobar_temperature = pmk_input%ref_isobar_temperature
                allocate(reference%isobar_volume( &
                    size(pmk_input%ref_isobar_volume)))
                reference%isobar_volume = 1.0e-3_dp*pmk_input%ref_isobar_volume
            end if
    
            reference%compare_density = pmk_input%compare_density
            if (reference%compare_density) then
                reference%density_weight = pmk_input%ref_density_weight
                reference%density = pmk_input%ref_density
                reference%density_temperature = pmk_input%ref_density_temperature
            end if
    
            reference%compare_phase_transition = pmk_input%compare_phase_transition
            if (reference%compare_phase_transition) then
                reference%phase_transition_weight = pmk_input%ref_phase_transition_weight
                reference%phase_transition = pmk_input%ref_phase_transition
            end if
        end subroutine qce_prepare
        !=================================================================================
        ! Initializes the conserved quantities (total number of particles and total mass
        ! of the system).
        subroutine initialize_conserved_quantities()
            integer:: i
    
            ! Particle number
            allocate(global_data%ntot(size(monomer)))
            global_data%ntot = global_data%monomer_amounts*avogadro
    
            ! Mass
            global_data%mtot = 0.0_dp
            do i = 1, size(monomer)
                global_data%mtot = global_data%mtot + &
                    global_data%monomer_amounts(i)*clusterset(monomer(i))%mass
            end do
            global_data%mtot = global_data%mtot / 1000.0_dp
    
            ! Excluded volume
            global_data%vexcl = 0.0_dp
            ! TODO: Generalize.
            select case (size(monomer))
            case (1)
                global_data%vexcl = global_data%monomer_amounts(1)*clusterset(monomer(1))%volume
            case (2)
                global_data%vexcl = &
                    global_data%bxv_pure(1)*global_data%monomer_amounts(1)*clusterset(monomer(1))%volume + &
                    global_data%bxv_pure(2)*global_data%monomer_amounts(2)*clusterset(monomer(2))%volume 
            end select
            global_data%vexcl = global_data%vexcl*avogadro*1.0e-30_dp
        end subroutine initialize_conserved_quantities
        !=================================================================================
        ! Initializes the degree of the population and mass polynomials.
        subroutine initialize_degree()
            integer:: i
            integer:: j
    
            global_data%degree = 0
            do i = 1, size(clusterset)
                do j = 1, size(monomer)
                    if (clusterset(i)%composition(j) > global_data%degree(j)) &
                        global_data%degree(j) = clusterset(i)%composition(j)
                end do
            end do
        end subroutine initialize_degree
        !=================================================================================
        ! Starts all QCE calculations. The main loop in here is OpenMP parallelized.
        subroutine qce_start()
            use auxiliary, only: progress_bar
            integer:: iamf
            integer:: ibxv
            integer:: itemp
            integer:: nr_isobars_computed, nr_isobars_total
    
            best_ib%error = huge(0.0)
            nr_isobars_computed = 0
            nr_isobars_total = global_data%amf%num*global_data%bxv%num

            !$OMP PARALLEL DEFAULT(none), &
            !$OMP& PRIVATE(ib, iamf, ibxv, itemp), &
            !$OMP& SHARED(global_data, clusterset, monomer, best_ib, reference, &
            !$OMP& nr_isobars_computed, nr_isobars_total)
    
            !$OMP DO COLLAPSE(2), SCHEDULE(GUIDED)
            do iamf = 1, global_data%amf%num
                do ibxv = 1, global_data%bxv%num
    
                    allocate(ib%temp(global_data%temp%num))
                    allocate(ib%vol(global_data%temp%num))
                    allocate(ib%gibbs(global_data%temp%num))
                    allocate(ib%converged(global_data%temp%num))
                    allocate(ib%solution(global_data%temp%num))
                    allocate(ib%lnq(global_data%temp%num, size(clusterset)))
                    allocate(ib%populations(global_data%temp%num, size(clusterset)))
    
                    ib%amf = global_data%amf%first + (iamf-1)*global_data%amf%delta
                    ib%amf = ib%amf/avogadro**2
                    ib%bxv = global_data%bxv%first + (ibxv-1)*global_data%bxv%delta
                    do itemp = 1, global_data%temp%num
                        ib%temp(itemp) = global_data%temp%first &
                            + (itemp-1)*global_data%temp%delta
                    end do
    
                    call qce_main(ib)
    
                    !$OMP CRITICAL
                    global_data%nconverged = global_data%nconverged + count(ib%converged)
                    if (ib%error < best_ib%error) best_ib = ib
                    !$OMP END CRITICAL
    
                    deallocate(ib%temp)
                    deallocate(ib%vol)
                    deallocate(ib%gibbs)
                    deallocate(ib%converged)
                    deallocate(ib%solution)
                    deallocate(ib%lnq)
                    deallocate(ib%populations)
    
                    !$OMP ATOMIC
                    nr_isobars_computed = nr_isobars_computed + 1
                    !$OMP END ATOMIC
    
#ifdef _OPENMP
                    if (omp_get_thread_num() == 0) then
                        call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar)
                    end if
#else
                        call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar)
#endif
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar, newline=.true.)
        end subroutine qce_start
        !=================================================================================
        ! Calculates a QCE isobar.
        subroutine qce_main(ib)
            type(isobar_t), intent(inout) :: ib
    
            real(dp):: v0
            real(dp):: vdamp
            real(dp):: vol
            real(dp):: gibbs
            real(dp):: error
            logical:: converged
            logical:: copy
            type(pf_t), dimension(size(clusterset)) :: lnq
            real(dp), dimension(size(clusterset)) :: populations
    
            integer:: itemp
    
            ! We perform two full QCE cycles for each temperature. In cycle One, we start
            ! from the ideal gas volume at the highest temperature. We decrease
            ! temperature and use the last converged volume as initial guess.
            ! In cycle Two, we start from a very small volume at the lowest
            ! temperature. We increase the temperature and use the last converged volume
            ! as initial guess.
            ! We choose the solution that led to smaller Gibbs free enthalpy.
    
            ! Cycle One
            vdamp = 1.0_dp - global_data%vdamp
            converged = .false.
            do itemp = size(ib%temp), 1, -1
                ! Use volume from previous temperature as initial guess, if available.
                ! Otherwise use ideal gas volume.
                if (.not. converged) then
                    v0 = ideal_gas_volume(ib%temp(itemp))
                else
                    v0 = vol
                end if
                ! Perform QCE iteration.
                call qce_iteration(ib%amf, ib%bxv, ib%temp(itemp), v0, vdamp, vol, &
                    gibbs, populations(:), lnq(:), 1, converged)
                ! Copy results.
                ib%gibbs(itemp) = gibbs
                ib%vol(itemp) = vol
                ib%lnq(itemp, :) = lnq(:)
                ib%converged(itemp) = converged
                ib%populations(itemp, :) = populations(:)
                ib%solution(itemp) = 1 ! Note below
                if(converged) ib%solution(itemp) = ib%solution(itemp) + 100 ! Note below.
            end do
    
            ! Cycle Two
            vdamp = 1.0_dp + global_data%vdamp
            converged = .false.
            do itemp = 1, size(ib%temp)
                ! Use volume from previous temperature as initial guess, if available.
                ! Otherwise use damped ideal gas volume.
                if (.not. converged) then
                    v0 = 1.0e-2_dp*ideal_gas_volume(ib%temp(itemp))
                else
                    !v0 = vol
                    ! Test:
                    v0 = min(vol, 1.0e-1*ideal_gas_volume(ib%temp(itemp)))
                end if
                ! Perform QCE iteration.
                call qce_iteration(ib%amf, ib%bxv, ib%temp(itemp), v0, vdamp, vol, &
                    gibbs, populations(:), lnq(:), 2, converged)
                ! Copy results, if necessary.
                if (converged) then
                    ib%solution(itemp) = ib%solution(itemp) + 10 ! Note below.
                    if (.not. ib%converged(itemp)) then
                        copy = .true.
                    else if (gibbs < ib%gibbs(itemp)) then
                        copy = .true.
                    else
                        copy = .false.
                    end if
                    if (copy) then
                        ib%gibbs(itemp) = gibbs
                        ib%vol(itemp) = vol
                        ib%lnq(itemp, :) = lnq(:)
                        ib%converged(itemp) = converged
                        ib%populations(itemp, :) = populations(:)
                        ib%solution(itemp) = ib%solution(itemp) + 1 ! Note below
                    end if
                end if
            end do
    
            ! Determine isobar quality, if necessary.
            ib%error = 0.0_dp
            if (reference%compare_density) then
                call compare_density(ib, reference%density_temperature, &
                    reference%density, error)
                ib%error = ib%error + reference%density_weight*error
            end if
            if (reference%compare_isobar) then
                call compare_isobar(ib, reference%isobar_temperature, &
                    reference%isobar_volume, error)
                ib%error = ib%error + reference%isobar_weight*error
            end if
            if (reference%compare_phase_transition) then
                call compare_phase_transition(ib, reference%phase_transition, error)
                ib%error = ib%error + reference%phase_transition_weight*error
            end if
    
            ! Note on ib%solution.
            ! This array indicates which solutions have converged and which solution was
            ! chosen. The first digit equals 1, if solution 1 has converged. The second
            ! digit equals 1, if solution 2 has converged. The third digits indicates
            ! which solution was chosen.
        end subroutine qce_main
        !=================================================================================
        ! Calculates the ideal gas volume at the given temperature.
        function ideal_gas_volume(temp)
            real(dp), intent(in) :: temp
            real(dp):: ideal_gas_volume
    
            ! The ideal gas volume.
            ideal_gas_volume = sum(global_data%ntot)*kb*temp/global_data%press
        end function ideal_gas_volume
        !=================================================================================
        ! Performs multiple QCE iterations and checks for convergence.
        subroutine qce_iteration(amf, bxv, temp, v0, vdamp, vol, gibbs, populations, &
            lnq, cyclus, converged)
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: v0
            real(dp), intent(in) :: vdamp
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: vol
            real(dp), intent(out) :: gibbs
            type(pf_t), dimension(size(clusterset)), intent(out) :: lnq
            real(dp), dimension(size(clusterset)), intent(out) :: populations
            integer, intent(in) :: cyclus
            logical, intent(out) :: converged
    
            integer:: iteration
            logical:: success
            real(dp):: vdamp_local
            real(dp):: old_vol
    
            ! Initialize the volume, populations, and Gibbs free enthalpy.
            vol = v0
            vdamp_local = vdamp
            gibbs = huge(0.0)
            call initialize_populations(populations)

            converged = .false.
            qce_loop: do iteration = 1, global_data%qce_iterations
                ! Calculate the cluster partition functions. In the first iteration, all
                ! of them need to be calculated. Later, we only need to update those,
                ! that depend on the volume.
                if (iteration == 1) then
                    call calculate_lnq(lnq, amf, bxv, temp, vol)
                else
                    call update_lnq(lnq, amf, bxv, temp, vol)
                end if
    
                ! Calculate new populations.
                call calculate_populations(populations, lnq, cyclus, success)
                if (.not. success) then
                    ! Without new populations, we can't proceed to get a new volume. In
                    ! order to proceed anyway, we will use a slightly damped initial
                    ! volume guess. The degree of damping depends on the number of times
                    ! that something in the QCE iteration has failed already.
                    vol = v0*vdamp_local
                    vdamp_local = vdamp_local*vdamp
                    call initialize_populations(populations)
                    cycle qce_loop
                end if
    
                ! Calculate new volume.
                old_vol = vol
                call calculate_volume(vol, vdamp, amf, bxv, temp, populations, success)
                if (.not. success) then
                    ! We couldn't get any physical volume and can't proceed to
                    ! recalculate the partition functions. In order to proceed anyway, we
                    ! will use a slightly damped initial volume guess. The degree of
                    ! damping depends on the number of times that something in the QCE
                    ! iteration has failed already.
                    vol = v0*vdamp_local
                    vdamp_local = vdamp_local*vdamp
                    call initialize_populations(populations)
                    cycle qce_loop
                end if
    
                ! Check for convergence.
                call check_convergence(gibbs, temp, vol, populations, lnq, converged)
                if (converged) exit qce_loop
            end do qce_loop
        end subroutine qce_iteration
        !=================================================================================
        ! Initializes populations. Assumes monomers, only.
        subroutine initialize_populations(populations)
            real(dp), dimension(size(clusterset)), intent(out) :: populations
    
            integer:: i
    
            ! Assume only monomers.
            populations = 0.0_dp
            do i = 1, size(monomer)
                populations(monomer(i)) = global_data%monomer_amounts(i)*avogadro
            end do
        end subroutine initialize_populations
        !=================================================================================
        ! Calculates populations, by solving the corresponding polynomials.
        ! Note that the coefficients of the polynomials in this subroutine are multiplied
        ! by the term (N_1^tot+N_2^tot)**(i+j) for numerical reasons. Thus the results of
        ! the Newton algorithm are populations divided by the same factor.
        subroutine calculate_populations(populations, lnq, cyclus, success)
            use polynomial
            real(dp), dimension(size(clusterset)), intent(out) :: populations
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            logical, intent(out) :: success
            integer, intent(in) :: cyclus
    
            integer:: iclust
            integer:: icomponent
            integer:: n
            integer:: m
            real(dp):: coeff
            real(dp):: sum_of_composition
            real(dp), dimension(size(monomer)) :: monomer_populations
            real(dp), dimension(size(clusterset)) :: pop_coeffs
            real(dp), dimension(size(clusterset)) :: mass_coeffs
            real(dp), dimension(:), allocatable :: coeffs1_pop
            real(dp), dimension(:, :), allocatable :: coeffs2_pop
            real(dp), dimension(:, :), allocatable :: coeffs2_mass
    
            ! Calculate coefficients for >> each cluster <<
            do iclust = 1, size(clusterset)
                associate(c => clusterset(iclust))
                    sum_of_composition = real(sum(c%composition), dp)
    
                    ! TEST
                    ! Calculate: ln[q(i)/(q(1)^i*q(2)^j)] + (i+j)*ln(N_1^tot+N_2^tot)
                    coeff = lnq(iclust)%qtot
                    do icomponent = 1, size(monomer)
                        coeff = coeff - &
                            c%composition(icomponent)*lnq(monomer(icomponent))%qtot
                    end do
                    coeff = coeff + sum_of_composition*log(sum(global_data%ntot))
    
                    ! Calculate coefficients of the population polynomial.
                    pop_coeffs(iclust) = exp(coeff)*sum_of_composition/sum(global_data%ntot)
                    ! Calculate coefficients of the mass polynomial.
                    mass_coeffs(iclust) = exp(coeff)*c%mass*amu/global_data%mtot
                end associate
            end do
    
            ! Calculate coefficients for >> each possible composition <<
            select case (size(monomer))
                case (1)
                    n = global_data%degree(1)
                    allocate(coeffs1_pop(0:n))
                    coeffs1_pop = 0.0_dp
                    coeffs1_pop(0) = -1.0_dp
                    do iclust = 1, size(clusterset)
                        associate(i => clusterset(iclust)%composition(1))
                            coeffs1_pop(i) = coeffs1_pop(i) + pop_coeffs(iclust)
                        end associate
                    end do
                case (2)
                    n = global_data%degree(1)
                    m = global_data%degree(2)
                    allocate(coeffs2_pop(0:n, 0:m))
                    allocate(coeffs2_mass(0:n, 0:m))
                    coeffs2_pop = 0.0_dp
                    coeffs2_pop(0, 0) = -1.0_dp
                    coeffs2_mass = 0.0_dp
                    coeffs2_mass(0, 0) = -1.0_dp
                    do iclust = 1, size(clusterset)
                        associate(i => clusterset(iclust)%composition(1), &
                            j => clusterset(iclust)%composition(2))
                            coeffs2_pop(i, j) = coeffs2_pop(i, j) + pop_coeffs(iclust)
                            coeffs2_mass(i, j) = coeffs2_mass(i, j) + mass_coeffs(iclust)
                        end associate
                    end do
            end select
    
            ! Initial guess of monomer populations.
            if (cyclus == 1) then
                ! Gas phase cycle.
                monomer_populations(:) = &
                    global_data%monomer_amounts(:)*avogadro/sum(global_data%ntot)
            else
                ! Liquid phase cycle.
                monomer_populations(:) = 1.0e-2_dp * &
                    global_data%monomer_amounts(:)*avogadro/sum(global_data%ntot)
            end if
            ! Solve the polynomials.
            select case (size(monomer))
                case (1)
                    call newton1(n, coeffs1_pop, monomer_populations(1), &
                        global_data%newton_iterations, success)
                    ! Check for unphysical solution.
                    if (monomer_populations(1) < 0.0_dp) success = .false.
                    deallocate(coeffs1_pop)
                case (2)
                    call newton2(n, m, coeffs2_pop, coeffs2_mass, monomer_populations, &
                        global_data%newton_iterations, success)
                    ! Check for unphysical solution.
                    if (any(monomer_populations < 0.0_dp)) success = .false.
                    deallocate(coeffs2_pop)
                    deallocate(coeffs2_mass)
            end select
    
            ! Calculate the remaining populations.
            if (success) then
                ! Get rid of the scaling factor N^tot and calculate the remaining
                ! populations.
                monomer_populations = monomer_populations*sum(global_data%ntot)
                call calculate_remaining_populations(populations, monomer_populations, &
                    lnq)
            end if
        end subroutine calculate_populations
        !=================================================================================
        ! Given the monomer populations, this calculates all other populations.
        subroutine calculate_remaining_populations(populations, monomer_populations, &
            lnq)
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            real(dp), dimension(size(monomer)), intent(in) :: monomer_populations
            real(dp), dimension(size(clusterset)), intent(out) :: populations
    
            integer:: i
            integer:: j
            real(dp):: tmp
    
            ! Assign monomer populations.
            do i = 1, size(monomer)
                populations(monomer(i)) = monomer_populations(i)
            end do
            ! Calculate the remaining populations.
            do i = 1, size(clusterset)
                associate(c => clusterset(i))
                    if (c%monomer) cycle
                    tmp = 0.0_dp
                    do j = 1, size(monomer)
                        tmp = tmp + real(c%composition(j), dp)* &
                            (log(monomer_populations(j)) - lnq(monomer(j))%qtot)
                    end do
                    tmp = tmp + lnq(i)%qtot
                    populations(i) = exp(tmp)
                end associate
            end do
        end subroutine calculate_remaining_populations
        !=================================================================================
        ! This subroutine calculates the new volume.
        subroutine calculate_volume(vol, vdamp, amf, bxv, temp, populations, success)
            use polynomial
            real(dp), intent(in) :: vdamp
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: vol
            real(dp), dimension(size(clusterset)), intent(in) :: populations
            logical, intent(out) :: success
    
            real(dp), dimension(0:3) :: coeffs
            complex(dp), dimension(3) :: roots
            logical, dimension(3) :: valid_roots

            real(dp), dimension(size(monomer)) :: amf_pure_local
            
            ! Choose old-style amf_pure, if not specified in input, i.e., if it is zero
            if (all(global_data%amf_pure <= 0.0_dp)) then
                amf_pure_local = amf
            else
                amf_pure_local = global_data%amf_pure
            end if

            ! TODO: Generalize
            select case (size(monomer))
            case(1)
                coeffs(0) = amf*bxv*global_data%vexcl*sum(global_data%ntot)**2
                coeffs(1) = -amf*sum(global_data%ntot)**2
            case(2)
                coeffs(0) = (global_data%ntot(1)**2*amf_pure_local(1) + &
                    2.0_dp*global_data%ntot(1)*global_data%ntot(2)*amf + &
                    global_data%ntot(2)**2*amf_pure_local(2))*bxv*global_data%vexcl
                coeffs(1) = -(global_data%ntot(1)**2*amf_pure_local(1) + &
                    2.0_dp*global_data%ntot(1)*global_data%ntot(2)*amf + &
                    global_data%ntot(2)**2*amf_pure_local(2))
            end select
    
            ! Calculate the coefficients.
            coeffs(2) = kb*temp*sum(populations) + &
                global_data%press*bxv*global_data%vexcl
            coeffs(3) = -global_data%press
    
            ! Solve the volume polynomial.
            call solve_polynomial3(coeffs, roots)
    
            ! Check for physical roots.
            where (abs(aimag(roots)) <= global_eps .and. &
                (real(roots, dp) - bxv*global_data%vexcl) >= global_eps)
                valid_roots = .true.
            else where
                valid_roots = .false.
            end where
    
            ! Use the smallest/largest valid volume, depending on the temperature cycle.
            if (any(valid_roots)) then
                if (vdamp <= 1.0_dp) then
                    vol = maxval(real(roots, dp), mask = valid_roots)
                else
                    vol = minval(real(roots, dp), mask = valid_roots)
                end if
                success = .true.
            else
                success = .false.
            end if
        end subroutine calculate_volume
        !=================================================================================
        ! This subroutine checks whether the QCE iteration has converged. It checks the
        ! deviation of the Gibbs free enthalpy from the previous run. The Gibbs free
        ! enthalpy ensures that both populations and volume have converged.
        subroutine check_convergence(gibbs, temp, vol, populations, lnq, &
            converged)
            use auxiliary, only: ln_factorial
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: temp
            real(dp), intent(inout) :: gibbs
            real(dp), dimension(size(clusterset)), intent(in) :: populations
            logical, intent(out) :: converged
            real(dp):: new_gibbs
            real(dp):: deviation
            integer:: iclust
    
            converged = .false.
            ! Calculate new Gibbs energy.
            new_gibbs = 0.0_dp
            do iclust = 1, size(clusterset)
                new_gibbs = new_gibbs - ln_factorial(populations(iclust)) + &
                    populations(iclust)*lnq(iclust)%qtot
            end do
            new_gibbs = -kb*temp*new_gibbs + global_data%press*vol
            deviation = new_gibbs/gibbs-1.0_dp
            gibbs = new_gibbs
    
            if (abs(deviation) <= global_data%max_deviation) converged = .true.
        end subroutine check_convergence
        !=================================================================================
        ! Compares to an experimental isobar.
        subroutine compare_isobar(ib, temp, vol, error)
            type(isobar_t), intent(in) :: ib
            real(dp), dimension(:), intent(in) :: temp
            real(dp), dimension(:), intent(in) :: vol
            real(dp), intent(out) :: error
    
            integer:: i
            real(dp):: vcalc
    
            error = 0.0_dp
            do i = 1, size(temp)
                vcalc = determine_volume(ib, temp(i))
                error = error + ((vcalc - vol(i))/vol(i))**2
            end do
            error = error / real(size(temp), dp)
        end subroutine compare_isobar
        !=================================================================================
        ! Compare to an experimental density at a given temperature.
        subroutine compare_density(ib, temp, density, error)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: density
            real(dp), intent(out) :: error
    
            real(dp):: vol
            real(dp):: vcalc
    
            vol = global_data%mtot/(1.0e3_dp*density)
            vcalc = determine_volume(ib, temp)
            error = ((vcalc - vol)/vol)**2
        end subroutine compare_density
        !=================================================================================
        ! Compares to an experimental phase transition.
        subroutine compare_phase_transition(ib, pt, error)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: pt
            real(dp), intent(out) :: error
    
            real(dp):: ptcalc
    
            ptcalc = determine_phase_transition(ib)
            error = ((ptcalc - pt)/pt)**2
        end subroutine compare_phase_transition
        !=================================================================================
        ! Determines the volume of a QCE isobar at a given temperature.
        ! Interpolates linearly.
        function determine_volume(ib, temp)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: temp
            real(dp):: determine_volume
    
            integer:: itemp
            real(dp):: x
    
            ! Narrow down temperature frame.
            search: do itemp = 1, size(ib%temp)
                if (ib%temp(itemp) >= temp) exit search
            end do search
    
            ! Temperature is between [itemp, itemp+1). Calculate interpolation factor.
            x = (temp - ib%temp(itemp))/(ib%temp(itemp+1) - ib%temp(itemp))
    
            ! Interpolate volume.
            determine_volume = (1.0_dp - x)*ib%vol(itemp) + x*ib%vol(itemp+1)
        end function determine_volume
        !=================================================================================
        ! Determines the phase transition of a QCE isobar, empirically by looking for the
        ! largest volume jump. This may not actually be a phase transition.
        function determine_phase_transition(ib)
            type(isobar_t), intent(in) :: ib
            real(dp):: determine_phase_transition
    
            integer:: itemp
            integer:: itemp_max
            real(dp):: dv
            real(dp):: dv_max
    
            dv_max = 0.0_dp
            itemp_max = 2
            do itemp = 2, size(ib%temp)
                if (ib%converged(itemp) .and. ib%converged(itemp-1)) then
                    dv = ib%vol(itemp) - ib%vol(itemp-1)
                    if (dv > dv_max) then
                        dv_max = dv
                        itemp_max = itemp
                    end if
                end if
            end do
    
            determine_phase_transition = &
                0.5_dp*(ib%temp(itemp_max) + ib%temp(itemp_max-1))
        end function determine_phase_transition
        !=================================================================================
        ! This subroutine prints info about the QCE calculation, calculates properties
        ! and writes output files.
        subroutine qce_finalize()
            ! Determine number of iterations and number of converged iterations both
            ! total and for the best isobar.
            write(*,'(4X,A,1X,G0,A,G0)') "Number of converged iterations:", &
                count(best_ib%converged), "/", global_data%temp%num
            if (reference%compare) then
                write(*,'(4X,A,1X,G0,A,G0)') &
                    "Number of converged iterations (total):", &
                    global_data%nconverged, "/", &
                    global_data%amf%num*global_data%bxv%num*global_data%temp%num
                write(*,*)
                ! TODO: Generalize.
                select case (size(monomer))
                case (1)
                    write(*,'(4X,A,G0.6,A,G0.6)') "Best isobar found for: amf = ", &
                        best_ib%amf*avogadro**2, " J*m^3/mol^2, bxv = ", &
                        best_ib%bxv
                case (2)
                    write(*,'(4X,A,G0.6,A,G0.6)') "Best isobar found for: amf_mix = ", &
                        best_ib%amf*avogadro**2, " J*m^3/mol^2, bxv_mix = ", &
                        best_ib%bxv
                end select
                write(*,'(4X,A,G0.6)') "Error: ", best_ib%error
                if (reference%compare_phase_transition) write(*,'(4X,A,G0.6,A)') &
                    "Calculated phase transition: ", &
                    determine_phase_transition(best_ib), " K"
                if (reference%compare_density) write(*,'(4X,A,G0.6,A)') &
                    "Calculated density: ", &
                    1.0e3_dp*global_data%mtot / &
                    (1.0e6_dp*determine_volume(best_ib, reference%density_temperature)), &
                    " g/cm^3"
            end if
            write(*,*)
    
            ! Perform Post-Processing.
            call post_processing(global_data%temp%num, size(clusterset), clusterset(:)%label, best_ib%temp, &
                best_ib%lnq, best_ib%populations, best_ib%vol, &
                best_ib%bxv*global_data%vexcl, best_ib%amf, global_data%press, &
                sum(global_data%ntot), best_ib%solution, best_ib%converged)
        end subroutine qce_finalize
        !=================================================================================
        ! This is the place, where all thermodynamic properties are caclulated using
        ! (hopefully) easy-to understand code. The idea is that this subroutine gets
        ! a couple of arrays containg parition functions, volume, populations, etc.
        ! and calculates everything there is to calculate from these in one place.
        ! Results are also written to file here.
        subroutine post_processing(ntemp, nclust, labels, temp, lnq_clust, pop, vol, vexcl, amf, &
            press, ntot, solution, converged)
            use auxiliary, only: ln_factorial, derivative
            use lengths
            use iso_varying_string
            use thermo
            integer, intent(in) :: ntemp
            integer, intent(in) :: nclust
            real(dp), intent(in) :: press
            real(dp), intent(in) :: ntot
            real(dp), intent(in) :: vexcl
            type(pf_t), dimension(ntemp, nclust) :: lnq_clust
            logical, dimension(ntemp), intent(in) :: converged
            integer, dimension(ntemp), intent(in) :: solution
            real(dp), dimension(ntemp), intent(in) :: vol
            real(dp), dimension(ntemp), intent(in) :: temp
            real(dp), intent(in) :: amf
            real(dp), dimension(ntemp, nclust), intent(in) :: pop
            type(varying_string), dimension(nclust), intent(in) :: labels

            integer:: iclust
            integer:: itemp
            integer:: myunit
            type(pf_t), dimension(ntemp, nclust) :: dlnq_clust, ddlnq_clust
            real(dp), dimension(ntemp) :: dvol
            real(dp), dimension(ntemp, nclust) :: pop_norm
            real(dp), dimension(ntemp, nclust) :: conc
            type(pf_t), dimension(nclust) :: d, dd
            character(fmt_len) :: fmtspec

            ! Calculate derivatives.
            do itemp = 1, ntemp
                call calculate_dlnq(d, amf, temp(itemp), vol(itemp))
                call calculate_ddlnq(dd, amf, temp(itemp), vol(itemp))
                dlnq_clust(itemp, :) = d
                ddlnq_clust(itemp, :) = dd
            end do
            dvol = derivative(vol, temp)

            ! Calculate and write thermodynamic quantities.
            call calculate_thermo(ntemp, nclust, pop, press, temp, vol, vexcl, &
                lnq_clust, dlnq_clust, ddlnq_clust, dvol, solution, converged)

            ! Normalize populations.
            do iclust = 1, nclust
                pop_norm(:, iclust) = &
                    pop(:, iclust)*real(sum(clusterset(iclust)%composition), dp)/ntot
            end do
    
            ! Write populations.
            write(fmtspec, '(A,G0,A)') '(A1,A12,', nclust, '(1X,A13))' 
            open(newunit = myunit, action = "write", status = "unknown", &
                file = "populations.dat")
            write(myunit, fmtspec) "#", "T/K", (char(labels(iclust)), iclust = 1, nclust)
            write(fmtspec, '(A,G0,A)') '(ES13.6,', nclust, '(1X,ES13.6))' 
            do itemp = 1, ntemp
                if (converged(itemp)) write(myunit, fmtspec) temp(itemp), &
                    (pop_norm(itemp, iclust), iclust = 1, nclust)
            end do
            close(myunit)
    
            ! Calculate concentrations.
            do iclust = 1, nclust
                conc(:, iclust) = pop(:, iclust)/(avogadro*vol*1.0e3_dp)
            end do
    
            ! Write concentrations.
            write(fmtspec, '(A,G0,A)') '(A1,A12,', nclust, '(1X,A13))' 
            open(newunit = myunit, action = "write", status = "unknown", &
                file = "concentrations.dat")
            write(myunit, fmtspec) "#", "T/K", (char(labels(iclust)), iclust = 1, nclust)
            write(fmtspec, '(A,G0,A)') '(ES13.6,', nclust, '(1X,ES13.6))' 
            do itemp = 1, ntemp
                if (converged(itemp)) write(myunit, fmtspec) temp(itemp), &
                    (conc(itemp, iclust), iclust = 1, nclust)
            end do
            close(myunit)
        end subroutine post_processing
        !=================================================================================
end module qce
!=========================================================================================
