module MOL_info
implicit none
type shells
    integer    :: angMoment
    integer    :: nGauss
    real(8),allocatable    :: exponents(:)
    real(8),allocatable    :: contrCoeff(:)
end type shells


type atom
    integer    :: charge
    real(8)    :: coor(3)
    character  :: base*30
    integer    :: nsShell
    integer    :: npShell
    integer    :: ndShell
    integer    :: nfShell
    integer    :: ngShell
    integer    :: nShell
    integer    :: nconts
    real(8)    :: atmForce(3)
    type(shells),allocatable :: shell(:)
end type atom
type(atom),allocatable :: atoms(:) 

!   Number of atoms  Number of contractions   Number of  contractions(spherical)
!         V                   V                  V 
integer Natoms,              Nconts ,            ncontssph

!===================
!     Martrix      +
!===================

!                      Fock_a  Fock_b     Overlap    Density :  core_hamilton     In UPtriangle format
!                       V         V          V         V           V
real(8),allocatable  :: Fa(:,:),  Fb(:,:), S(:,:),   P(:,:),   Hcore(:,:)    


!                        Density_a   Density_b
!                        V           V
real(8),allocatable  ::  Pa(:,:),    Pb(:,:)

!                       2-electron integral                    
!                          V                         
real(8),allocatable  :: TWOEI(:)


!                     Coeffient_a   Coeffient_b
!                       V           V
real(8),allocatable  :: C_a(:,:),   C_b(:,:)

!        Energy level   alpha      beta 
!                        V         V
real(8),allocatable  :: eLev_a(:),eLev_b(:)


!                  Transform Matrix: othogonal
!                       V
real(8),allocatable  :: X(:,:)


!Basis Normalization
!
real(8),allocatable  :: NorVEC(:)
real(8),allocatable  :: NorVECsph(:)
!===================
!    Properties    +
!===================
integer              :: Charge        ! Total charge
integer              :: Multi         ! Multiplicity
real(8),allocatable  :: MulCharge(:)  ! Mulliken Charge
real(8),allocatable  :: Force(:,:)    ! Force
real(8)              :: E,E_rep       ! Energy and core repulsion energy

integer              :: n_alpha       ! alpha electrons
integer              :: n_beta        ! beta electrons
integer              :: coreChg       ! total core charges

integer,allocatable  :: linkMat(:,:)  ! Linkage matrix

!===================
!   Functional     +
!===================
character Functional*30

!===================
!   integral info  +
!===================
integer,allocatable :: atm(:,:)
integer,allocatable :: bas(:,:)
real(8),allocatable :: env(:)
integer             :: nBases
end module MOL_info


module GRID_info
!===================
!   Grid           +
!===================
implicit none
type Grid
    real(8)    :: coor(3)
    real(8)    :: weight
    real(8)    :: rho_a,rho_b
    real(8)    :: sigma_aa,sigma_ab,sigma_bb
    real(8),allocatable ::val0(:)
    real(8),allocatable ::val1(:,:)
    real(8),allocatable :: TempD(:)
end type Grid

integer :: ngrids

type(Grid),allocatable :: Grids(:) 
end module GRID_info



