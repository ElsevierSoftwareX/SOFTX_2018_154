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
! This module implements algorithms for solving polynomials.
!=========================================================================================
module polynomial
    use kinds
    use cluster
    use constants
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: newton1, newton2
    public :: solve_polynomial3
    !=====================================================================================
    ! Convergence criterion in the Newton algorithm.
    real(dp), parameter :: newton_convergence = global_eps
    !=====================================================================================
    contains
        !=================================================================================
        ! The Newton-Raphson algorithm for the solution of a 1D polynomial.
        subroutine newton1(n, c, x, iterations, success)
            integer, intent(in) :: n
            real(dp), dimension(n), intent(in) :: c
            real(dp), intent(inout) :: x
            integer, intent(in) :: iterations
            logical, intent(out) :: success

            integer :: i
            real(dp) :: p, pdiff, dx, xnew

            newton_loop: do i = 1, iterations
                call horner1(n, c, x, p, pdiff)

                ! Check for convergence.
                if (abs(p) <= newton_convergence) then
                    success = .true.
                    return
                end if

                ! Determine Newton step.
                dx = -p/pdiff

                ! Half step size, if this would move x out of physical bounds.
                xnew = x + dx
                do while (xnew < 0.0_dp .or. xnew > 1.0_dp)
                    dx = 0.5 * dx
                    xnew = x + dx
                end do

                ! Move solution closer to zero.
                x = xnew
            end do newton_loop
            success = .false.
        end subroutine newton1
        !=================================================================================
        ! The Newton-Raphson algorithm for the simultaneous solution of two 2D 
        ! polynomials.
        subroutine newton2(n, m, c1, c2, x, iterations, success)
            integer, intent(in) :: n, m
            real(dp), dimension(n, m), intent(in) :: c1, c2
            real(dp), dimension(2), intent(inout) :: x
            integer, intent(in) :: iterations
            logical, intent(out) :: success

            integer :: i, j
            real(dp), dimension(2) :: p, dx, xnew
            real(dp), dimension(2, 2) :: pdiff
            real(dp) :: det

            newton_loop: do i = 1, iterations
                call horner2(n, m, c1, x, p(1), pdiff(:, 1))
                call horner2(n, m, c2, x, p(2), pdiff(:, 2))
                pdiff = transpose(pdiff)
                
                ! Check for convergence.
                if (all(abs(p) <= newton_convergence)) then
                    success = .true.
                    return
                end if

                ! Determine Newton step.
                det = pdiff(1,1)*pdiff(2,2) - pdiff(1,2)*pdiff(2,1)
                dx(1) = p(2)*pdiff(1,2) - p(1)*pdiff(2,2)
                dx(2) = p(1)*pdiff(2,1) - p(2)*pdiff(1,1)
                dx = dx/det

                ! Half step size, if this would move x out of physical bounds.
                do j = 1, 2
                    xnew(j) = x(j) + dx(j)
                    do while (xnew(j) < 0.0_dp .or. xnew(j) > 1.0_dp)
                        dx = 0.5 * dx
                        xnew(j) = x(j) + dx(j)
                    end do
                end do

                ! Move solution closer to zero.
                x = xnew
            end do newton_loop
            success = .false.
        end subroutine newton2
        !=================================================================================
        ! Solves a third order polynomial.
        subroutine solve_polynomial3(coeffs, roots)
            real(dp), dimension(0:3), intent(in) :: coeffs
            complex(dp), dimension(3), intent(out) :: roots

            complex(dp), dimension(6) :: cof
            complex(dp), parameter :: im = (0.0_dp,1.0_dp)
            complex(dp), parameter :: one_third = cmplx(1.0_dp/3.0_dp, 0.0_dp, dp)
            real(dp), parameter :: two_power_one_third = 2.0_dp**(1.0_dp/3.0_dp)
            real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
            
            roots = (0.0_dp, 0.0_dp)
            cof = (0.0_dp, 0.0_dp)
            
            cof(1) = cmplx(-coeffs(2)/3.0_dp/coeffs(3), 0.0_dp, kind = dp)
            cof(2) = cmplx(-coeffs(2)**2 + 3.0_dp*coeffs(3)*coeffs(1), 0.0_dp, kind = dp)
            cof(3) = cmplx(-2.0_dp*coeffs(2)**3 + 9.0_dp*coeffs(3)*coeffs(2)*coeffs(1) &
                -27.0_dp*coeffs(0)*coeffs(3)**2, 0.0_dp, kind = dp)
            cof(4) = 4.0_dp*cof(2)**3 + cof(3)**2
            cof(5) = cmplx(3.0_dp*coeffs(3), 0.0_dp, kind = dp)
            cof(6) = (cof(3) + sqrt(cof(4)))**one_third
            
            roots(1) = cof(1) - two_power_one_third*cof(2)/(cof(5)*cof(6)) + &
                cof(6)/two_power_one_third/cof(5)
            roots(2) = cof(1) + &
                (1.0_dp + im*sqrt3)*cof(2)/(two_power_one_third**2*cof(5)*cof(6)) - &
                (1.0_dp - im*sqrt3)*cof(6)/(2.0_dp*cof(5)*two_power_one_third)
            roots(3) = cof(1) + &
                (1.0_dp - im*sqrt3)*cof(2)/(two_power_one_third**2*cof(5)*cof(6)) - &
                (1.0_dp + im*sqrt3)*cof(6)/(2.0_dp*cof(5)*two_power_one_third)
        end subroutine solve_polynomial3
        !=================================================================================
        ! Uses Horner's scheme to evaluate a 1D polynomial of degree n and its derivative
        ! at a point x, given the coefficients c0, ..., cn.
        subroutine horner1(n, c, x, p, pdiff)
            integer, intent(in) :: n
            real(dp), dimension(0:n), intent(in) :: c
            real(dp), intent(in) :: x
            real(dp), intent(out) :: p, pdiff

            integer :: i

            p = c(n)
            pdiff = 0.0_dp
            do i = n-1, 0, -1
                pdiff = x*pdiff + p
                p = x*p + c(i)
            end do
        end subroutine horner1
        !=================================================================================
        ! Uses Horner's scheme to evaluate a 2D polynomial of degree (n,m) and its
        ! derivative at a point x, given the coefficients c00, ..., cnm.
        subroutine horner2(n, m, c, x, p, pdiff)
            integer, intent(in) :: n, m
            real(dp), dimension(0:n, 0:m), intent(in) :: c
            real(dp), dimension(2), intent(in) :: x
            real(dp), intent(out) :: p
            real(dp), dimension(2), intent(out) :: pdiff

            integer :: j
            real(dp) :: b, bdiff

            call horner1(n, c(:, m), x(1), b, bdiff)
            p = b
            pdiff(1) = bdiff
            pdiff(2) = 0.0_dp
            do j = m-1, 0, -1
                call horner1(n, c(:, j), x(1), b, bdiff)
                pdiff(1) = x(2)*pdiff(1) + bdiff
                pdiff(2) = x(2)*pdiff(2) + p
                p = x(2)*p + b
            end do
        end subroutine horner2
        !=================================================================================
end module polynomial
!=========================================================================================
