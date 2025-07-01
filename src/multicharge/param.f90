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

module multicharge_param
   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type
   use mctc_io_convert, only: autoaa
   use mctc_data, only: get_covalent_rad, get_pauling_en, get_vdw_rad
   use multicharge_model, only: mchrg_model_type, &
      & new_eeq_model, eeq_model, new_eeqbc_model, eeqbc_model, &
      & new_eeqbceps_model, eeqbceps_model
   use multicharge_param_eeq2019, only: get_eeq_chi, get_eeq_eta, &
      & get_eeq_rad, get_eeq_kcnchi
   use multicharge_param_eeqbc2024, only: get_eeqbc_chi, get_eeqbc_eta, &
      & get_eeqbc_rad, get_eeqbc_kcnchi, get_eeqbc_kqchi, get_eeqbc_kqeta, &
      & get_eeqbc_cap, get_eeqbc_cov_radii, get_eeqbc_avg_cn
   use multicharge_param_eeqbceps2025, only: get_eeqbceps_chi, get_eeqbceps_eta, &
      & get_eeqbceps_rad, get_eeqbceps_kcnchi, get_eeqbceps_kqchi, &
      & get_eeqbceps_kqeta, get_eeqbceps_cap, get_eeqbceps_cov_radii, &
      & get_eeqbceps_avg_cn, get_eeqbceps_born_radii
   implicit none
   private

   public :: new_eeq2019_model, new_eeqbc2024_model, new_eeqbceps2025_model, mchargeModel

   !> Possible charge models enumerator
   type :: TMchargeModelEnum
      !> Classic electronegativity equilibration model
      integer :: eeq2019 = 1
      !> Bond-capacitor corrected electronegativity equilibration model
      integer :: eeqbc2024 = 2
      !> Bond-capacitor corrected electronegativity equilibration model with implicit solvation
      integer :: eeqbceps2025 = 3
   end type TMchargeModelEnum

   !> Actual charge model enumerator
   type(TMchargeModelEnum), parameter :: mchargeModel = TMchargeModelEnum()

contains

   subroutine new_eeq2019_model(mol, model, error)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Electronegativity equilibration model
      class(mchrg_model_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: cutoff = 25.0_wp
      real(wp), parameter :: cn_exp = 7.5_wp
      real(wp), parameter :: cn_max = 8.0_wp

      real(wp), allocatable :: chi(:), eta(:), kcnchi(:), rad(:), rcov(:)
      type(eeq_model), allocatable :: eeq

      chi = get_eeq_chi(mol%num)
      eta = get_eeq_eta(mol%num)
      kcnchi = get_eeq_kcnchi(mol%num)
      rad = get_eeq_rad(mol%num)
      rcov = get_covalent_rad(mol%num)

      allocate (eeq)
      call new_eeq_model(eeq, mol=mol, error=error, chi=chi, &
         & rad=rad, eta=eta, kcnchi=kcnchi, cutoff=cutoff, &
         & cn_exp=cn_exp, rcov=rcov, cn_max=cn_max)
      call move_alloc(eeq, model)

   end subroutine new_eeq2019_model

   subroutine new_eeqbc2024_model(mol, model, error)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Electronegativity equilibration model
      class(mchrg_model_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: kcnrad = 0.14_wp
      real(wp), parameter :: kbc = 0.60_wp
      real(wp), parameter :: cutoff = 25.0_wp
      real(wp), parameter :: cn_exp = 2.0_wp
      real(wp), parameter :: norm_exp = 0.75_wp

      real(wp), allocatable :: chi(:), eta(:), rad(:), kcnchi(:), &
         & kqchi(:), kqeta(:), cap(:), rcov(:), avg_cn(:), en(:), &
         & rvdw(:, :)
      type(eeqbc_model), allocatable :: eeqbc

      chi = get_eeqbc_chi(mol%num)
      eta = get_eeqbc_eta(mol%num)
      rad = get_eeqbc_rad(mol%num)
      kcnchi = get_eeqbc_kcnchi(mol%num)
      kqchi = get_eeqbc_kqchi(mol%num)
      kqeta = get_eeqbc_kqeta(mol%num)
      cap = get_eeqbc_cap(mol%num)
      rcov = get_eeqbc_cov_radii(mol%num)
      avg_cn = get_eeqbc_avg_cn(mol%num)
      ! Electronegativities normalized to Fluorine
      ! with actinides (Th-Lr) set to average of 1.30
      en = get_pauling_en(mol%num)
      en = merge(en, 1.30_wp, mol%num < 90)
      en = merge(0.80_wp, en, mol%num == 87)
      en = merge(1.00_wp, en, mol%num == 89)
      en = merge(1.10_wp, en, mol%num == 90 .or. mol%num == 91 &
         &.or. mol%num == 92 .or. mol%num == 95)
      en = merge(1.20_wp, en, mol%num == 93 .or. mol%num == 94 &
         &.or. mol%num == 97 .or. mol%num == 103)
      en = en/3.98_wp
      rvdw = get_vdw_rad(spread(mol%num(mol%id), 2, mol%nat), &
         & spread(mol%num(mol%id), 1, mol%nat)) * autoaa

      allocate (eeqbc)
      call new_eeqbc_model(eeqbc, mol=mol, error=error, chi=chi, & 
         & rad=rad, eta=eta, kcnchi=kcnchi, kqchi=kqchi, kqeta=kqeta, &
         & kcnrad=kcnrad, cap=cap, avg_cn=avg_cn, kbc=kbc, &
         & cutoff=cutoff, cn_exp=cn_exp, rcov=rcov, en=en, &
         & norm_exp=norm_exp, rvdw=rvdw)
      call move_alloc(eeqbc, model)

   end subroutine new_eeqbc2024_model

   subroutine new_eeqbceps2025_model(mol, model, error, epsilon, paramspath)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Electronegativity equilibration model
      class(mchrg_model_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Epsilon for the implicit Born model
      real(wp), optional, intent(in) :: epsilon
      !> Path to the parameters file
      character(len=*), optional, intent(in) :: paramspath

      real(wp), parameter :: kcnrad = 0.14_wp
      real(wp), parameter :: kbc = 0.60_wp
      real(wp), parameter :: cutoff = 25.0_wp
      real(wp), parameter :: cn_exp = 2.0_wp
      real(wp), parameter :: norm_exp = 0.75_wp

      real(wp), allocatable :: chi(:), eta(:), rad(:), kcnchi(:), &
         & kqchi(:), kqeta(:), cap(:), rcov(:), avg_cn(:), en(:), &
         & rvdw(:, :), radii(:), eps

      type(eeqbceps_model), allocatable :: eeqbceps

      if (present(paramspath)) then
         chi = read_param_file(paramspath // 'chi.txt', mol%num)
         eta = read_param_file(paramspath // 'eta.txt', mol%num)
         rad = read_param_file(paramspath // 'rad.txt', mol%num)
         kcnchi = read_param_file(paramspath // 'kcnchi.txt', mol%num)
         kqchi = read_param_file(paramspath // 'kqchi.txt', mol%num)
         kqeta = read_param_file(paramspath // 'kqeta.txt', mol%num)
         cap = read_param_file(paramspath // 'cap.txt', mol%num)
         rcov = read_param_file(paramspath // 'rcov.txt', mol%num)
         radii = read_param_file(paramspath // 'radii.txt', mol%num)
      else
         chi = get_eeqbceps_chi(mol%num)
         eta = get_eeqbceps_eta(mol%num)
         rad = get_eeqbceps_rad(mol%num)
         kcnchi = get_eeqbceps_kcnchi(mol%num)
         kqchi = get_eeqbceps_kqchi(mol%num)
         kqeta = get_eeqbceps_kqeta(mol%num)
         cap = get_eeqbceps_cap(mol%num)
         rcov = get_eeqbceps_cov_radii(mol%num)
         radii = get_eeqbceps_born_radii(mol%num)
      end if

      avg_cn = get_eeqbceps_avg_cn(mol%num)
      
      if (present(epsilon)) then
         ! Use the provided epsilon value
         eps = epsilon
      else
         ! Use the default epsilon value
         eps = 1.0_wp
      end if

      ! Electronegativities normalized to Fluorine
      ! with actinides (Th-Lr) set to average of 1.30
      en = get_pauling_en(mol%num)
      en = merge(en, 1.30_wp, mol%num < 90)
      en = merge(0.80_wp, en, mol%num == 87)
      en = merge(1.00_wp, en, mol%num == 89)
      en = merge(1.10_wp, en, mol%num == 90 .or. mol%num == 91 &
         &.or. mol%num == 92 .or. mol%num == 95)
      en = merge(1.20_wp, en, mol%num == 93 .or. mol%num == 94 &
         &.or. mol%num == 97 .or. mol%num == 103)
      en = en/3.98_wp
      rvdw = get_vdw_rad(spread(mol%num(mol%id), 2, mol%nat), &
         & spread(mol%num(mol%id), 1, mol%nat)) * autoaa

      allocate (eeqbceps)
      call new_eeqbceps_model(eeqbceps, mol=mol, error=error, chi=chi, & 
         & rad=rad, eta=eta, kcnchi=kcnchi, kqchi=kqchi, kqeta=kqeta, &
         & kcnrad=kcnrad, cap=cap, avg_cn=avg_cn, kbc=kbc, &
         & cutoff=cutoff, cn_exp=cn_exp, rcov=rcov, en=en, &
         & norm_exp=norm_exp, rvdw=rvdw, &
         & radii=radii, epsilon=eps)
      call move_alloc(eeqbceps, model)

   end subroutine new_eeqbceps2025_model


   function read_param_file(filename, num) result(param)
      !> Read a parameter file and return the parameters
      character(len=*), intent(in) :: filename
      !> Atomic numbers of the elements
      integer, intent(in) :: num(:)
      !> Parameters read from the file
      real(wp), allocatable :: param(:), raw_param(:)


      integer :: i, n, io
      character(len=256) :: line
      open(unit=10, file=filename, action='read')
      read(10, '(A)') line
      allocate(raw_param(103))
      do i = 1, 103
         read(10, '(F20.10)', iostat=io) raw_param(i)
         if (i /= 103 .and. io /= 0) then
            close(10)
            error stop 'Error reading parameter file: '//trim(filename)
         end if
      end do
      close(10)

      n = size(num)
      allocate(param(n))
      do i = 1, n
         if (num(i) < 1 .or. num(i) > 103) then
            error stop 'Invalid atomic number:'
         end if
         param(i) = raw_param(num(i))
      end do

   end function read_param_file


end module multicharge_param
