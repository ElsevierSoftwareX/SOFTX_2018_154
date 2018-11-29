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
! This module provides atomic data.
module atomic_data
    use kinds
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: periodic_table
    !=====================================================================================
    ! Data type storing atomic data.
    type :: element_t
        integer :: atomic_number
        character(2) :: atomic_symbol
        real(dp) :: atomic_mass ! in amu
    end type element_t
    !=====================================================================================
    ! Class that provides access to the periodic table of elements.
    type :: periodic_table_t
        integer :: n
        type(element_t), dimension(:), allocatable :: elements
    contains
        procedure, public :: init => periodic_table_init
        procedure, public :: mass => periodic_table_mass
    end type periodic_table_t
    !=====================================================================================
    ! Instantiation of the periodic_table_t type.
    type(periodic_table_t):: periodic_table
    !=====================================================================================
    contains
        !=================================================================================
        ! Initializes the periodic table.
        subroutine periodic_table_init(table)
            use constants
            class(periodic_table_t):: table

            table%n = 18
            allocate(table%elements(table%n))

            table%elements( 1) = element_t( 1, "H ",  1.008)
            table%elements( 2) = element_t( 2, "He",  4.003)
            table%elements( 3) = element_t( 3, "Li",  6.94 )
            table%elements( 4) = element_t( 4, "Be",  9.012)
            table%elements( 5) = element_t( 5, "B ", 10.81 )
            table%elements( 6) = element_t( 6, "C ", 12.01 )
            table%elements( 7) = element_t( 7, "N ", 14.01 )
            table%elements( 8) = element_t( 8, "O ", 16.00 )
            table%elements( 9) = element_t( 9, "F ", 19.00 )
            table%elements(10) = element_t(10, "Ne", 20.18 )
            table%elements(11) = element_t(11, "Na", 22.99 )
            table%elements(12) = element_t(12, "Mg", 24.31 )
            table%elements(13) = element_t(13, "Al", 26.98 )
            table%elements(14) = element_t(14, "Si", 28.09 )
            table%elements(15) = element_t(15, "P ", 30.97 )
            table%elements(16) = element_t(16, "S ", 32.06 )
            table%elements(17) = element_t(17, "Cl", 35.45 )
            table%elements(18) = element_t(18, "Ar", 39.95 )
        end subroutine periodic_table_init
        !=================================================================================
        ! Returns the mass of an element.
        function periodic_table_mass(table, atomic_symbol) result(mass)
            use error
            class(periodic_table_t), intent(in) :: table
            character(2), intent(in) :: atomic_symbol
            real(dp):: mass

            integer:: i
            do i = 1, table%n
                if (table%elements(i)%atomic_symbol == atomic_symbol) then
                    mass = table%elements(i)%atomic_mass
                    return
                end if
            end do

            mass = 0.0_dp
            call pmk_error("unknown element '" // atomic_symbol // "'")
        end function periodic_table_mass
        !=================================================================================
end module atomic_data
!=========================================================================================
