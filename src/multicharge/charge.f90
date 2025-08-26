! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/charge.f90
!> Contains functions to calculate the partial charges with and without
!> separate charge model setup

!> Interface to the charge models
module multicharge_charge
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_cutoff, only : get_lattice_points
   use multicharge_model, only : mchrg_model_type, eeqbceps_model
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model, new_eeqbceps2025_model
   use multicharge_param_eeqbceps2025, only: get_eeqbceps_rad, get_eeqbceps_avg_cn
   implicit none
   private

   public :: get_charges, get_eeq_charges, get_eeqbc_charges, get_eeqbceps_charges


contains


!> Classical electronegativity equilibration charges
subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL)

   !> Multicharge model
   class(mchrg_model_type), intent(in) :: mchrg_model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   logical :: grad
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: trans(:, :)

   grad = present(dqdr) .and. present(dqdL)

   allocate(cn(mol%nat), qloc(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      allocate (dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, mchrg_model%ncoord%cutoff, trans)
   call mchrg_model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call mchrg_model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call mchrg_model%solve(mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & qvec=qvec, dqdr=dqdr, dqdL=dqdL)

end subroutine get_charges


!> Obtain charges from electronegativity equilibration model
subroutine get_eeq_charges(mol, error, qvec, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   class(mchrg_model_type), allocatable :: eeq_model

   call new_eeq2019_model(mol, eeq_model, error)

   call get_charges(eeq_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeq_charges


!> Obtain charges from bond capacity electronegativity equilibration model
subroutine get_eeqbc_charges(mol, error, qvec, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   class(mchrg_model_type), allocatable :: eeqbc_model

   call new_eeqbc2025_model(mol, eeqbc_model, error)

   call get_charges(eeqbc_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeqbc_charges


!> Obtain charges from the epsilon dependent electronegativity equilibration model
subroutine get_eeqbceps_charges(mol, epsilon, error, qvec, aI, bornradscal, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dielectric permittivity of the medium
   real(wp), intent(in) :: epsilon

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   real(wp), intent(out), allocatable, optional :: ai(:)

   real(wp), intent(in), optional :: bornradscal

   class(mchrg_model_type), allocatable :: eeqbceps_model

  real(wp), allocatable :: cn(:)
  real(wp), allocatable :: trans(:,:)
  
  ! locals for optional ai computation
  real(wp), allocatable :: rad(:), avg_cn(:)
  real(wp) :: norm_cn, radi
  real(wp), parameter :: kcnrad = 0.14_wp
  real(wp), parameter :: norm_exp = 0.75_wp
  integer :: iat, izp

   call new_eeqbceps2025_model(mol, eeqbceps_model, error, epsilon=epsilon, &
      bornradscal=bornradscal)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, eeqbceps_model%ncoord%cutoff, trans)
   call eeqbceps_model%ncoord%get_coordination_number(mol, trans, cn)

   if (present(ai)) then
      ! Compute ai the same way as eeqbceps_model%getradiiding would do
      rad = get_eeqbceps_rad(mol%num)
      avg_cn = get_eeqbceps_avg_cn(mol%num)
      allocate(ai(mol%nat))
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / avg_cn(izp)**norm_exp
         radi = rad(izp) * (1.0_wp - kcnrad*cn(iat)*norm_cn)
         ai(iat) = radi
      end do
   end if

   call get_charges(eeqbceps_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeqbceps_charges

end module multicharge_charge
