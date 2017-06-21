

      MODULE usr


        Use param
        Use param1
        Use physprop
        USE run
        USE ic, only: IC_X_s
        USE compar, only: MyPE,PE_IO
        USE rxn_com

        use rxns, only: NO_OF_RXNS
        use rxns, only: Reaction


!
!       Declare the user-defined namelist variables (usrnlst.inc) in this module.
!       Also Include user-defined variables in this module.  To access the
!       variables from a subroutine add the statement "Use usr".  If allocatable
!       arrays are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!
!   Author: David Huckaby - huckaby@netl.doe.gov                          
!   Date: 6-Apr-12
!
!   Module calculate the reaction rates of NETL amine sorbent 196C and 32D.
!   Reaction kinetics are based on the work of Lee, Mebane and Fauth (32D) and 
!   Membane Lee and Fauth (196C)

!   Modified by: Jeff Dietiker



    implicit none

!   a dummy variable listed in usrnlst.inc
    DOUBLE PRECISION DUMMY_DP

    logical :: use_conc


    type unit_type
        real(kind=8) :: m3
        real(kind=8) :: mol
        real(kind=8) :: Pa
        real(kind=8) :: J
        character(len=5) :: name
    end type

    type (unit_type) :: unit_MKS, unit_SI, unit_cgs




    type reaction_type
        real(kind=8) :: P_ref
        real(kind=8) :: nv
        real(kind=8) :: Cbar
        real(kind=8) :: E 
        real(kind=8) :: dH    
        real(kind=8) :: dS
        real(kind=8) :: zeta
        real(kind=8) :: A
        real(kind=8) :: B
        real(kind=8) :: dM          ! net mass transfered to the gas from the solid in kg/kmol
        real(kind=8) :: m
    end type


    integer, parameter :: nReac = 3

! WARNING: The following indices must match those defined inm species.inc
    integer, parameter :: i_N2 = 1, i_CO2 = 2, i_H2O_g = 3
    integer, parameter :: i_R2NH = 2, i_R2NCO2 = 3, i_R2NH2 = 4, i_HCO3 = 5, i_H2O_s = 6

    real(kind=8), parameter :: Rgas_MKS = 8.314471e3, Rgas_SI = 8.314471, Rgas_cgs = 8.314471e7, Rgas_cal = 1.9858775
    real(kind=8) :: Rgas_energy, Rgas_press
    real(kind=8), parameter :: X_small = 1e-7, X_zero = 1e-20
    !real(kind=8), parameter :: unit_m3 = 1e6, unit_mol = 1., unit_Pa = 10., unit_J = 1e7    

    type (reaction_type), dimension(nReac), target :: reac   

    INTEGER :: RR

    type (unit_type) :: unit

    DOUBLE PRECISION, DIMENSION(nReac) :: REAC_E, REAC_DH, REAC_DS,REAC_LOGZETA

    contains

    subroutine unit_init

        ! converts unit to base-unit system 
        unit_cgs%m3  = 1e6
        unit_cgs%mol = 1.
        unit_cgs%Pa  = 10.0
        unit_cgs%J   = 1e7  
        unit_cgs%name = 'cgs'

        unit_MKS%m3  = 1.
        unit_MKS%mol = 1e-3
        unit_MKS%Pa  = 1.
        unit_MKS%J   = 1.
        unit_MKS%name = 'MKS'

        unit_SI%m3   = 1.
        unit_SI%mol  = 1.
        unit_SI%Pa   = 1.
        unit_SI%J    = 1.
        unit_SI%name = 'SI'

    end subroutine




    subroutine reaction_init(unitsystem_name)

        implicit none


        character(len=*), intent(in) :: unitsystem_name

        type(reaction_type), pointer :: r => NULL()
        integer :: iSpec, iReac
        real(kind=8), parameter :: tol = 1e-7
        real(kind=8) :: P0,nv0

        INCLUDE 'species.inc' 

        call unit_init()

        ! what MFIX calls SI is MKS
        if (unitsystem_name .eq. 'MFIX_SI') then
            Rgas_energy = Rgas_MKS
            Rgas_press  = Rgas_MKS
            unit = unit_MKS
        else if (unitsystem_name .eq. 'MFIX_CGS') then
            Rgas_energy = Rgas_cal   ! mfix uses cal for heat/energy
            Rgas_press  = Rgas_cgs
            unit = unit_cgs
        else
            write(*,*) 'Unit system requested = ',unitsystem_name, ' provided = cgs'
            Rgas_energy = Rgas_cal
            unit = unit_cgs
        end if

!       write(*,*) 'reaction init - start'



         ! nv is computed automatically from inputs in mfix.dat
         ! nv = (amine_weight_fraction x solids density) / Molecular weight of amine
         ! factor of 1000.0 used for unit conversion (mol/kg)
         IF(RUN_TYPE=='NEW') THEN
            nv0 = IC_X_s(1,1,2) * RO_s0(1)/MW_s(1,2) * 1000.0
            if(nv0 .eq. 0) then
                write(*,*) 'IC#1 does not contain solid phase 1, check IC#2'
                nv0 = IC_X_s(2,1,2) * RO_s0(1)/MW_s(1,2) * 1000.0
                if(nv0 .eq. 0) then
                        write(*,*) 'IC#2 does not contain solid phase 1, check IC#3'
                        nv0 = IC_X_s(3,1,2) * RO_s0(1)/MW_s(1,2) * 1000.0
                end if
            end if
         ELSE              ! Read value of nv0 from file when doing a restart
            open (unit=888,file='nv0.dat')
            read(888,*)nv0
            close(888)
         ENDIF

         P0 = 97189.3

         ! 1. Dry CO2 absorption
         ! 2R2NH + CO2 <--> R2NH2 + R2NCO2
         ! This correspond to equations 1 (forward) and 2 (reverse)
         ! equations in mfix.dat

         r => reac(1)

         r%P_ref = P0 * unit%Pa                  ! dynes/cc        
         r%nv   = nv0 * unit%mol/unit%m3         ! mol-sites/m3
         r%Cbar = r%nv                           ! mol/cc - 1 mol-R2NH3/mol-site

         r%E    = REAC_E(1)/Rgas_SI              ! K  
         r%dH   = REAC_DH(1)/Rgas_SI             ! K
         r%dS   = REAC_DS(1)/Rgas_SI             
         r%zeta = 10**REAC_LOGZETA(1)/unit%Pa   

         r%A = r%zeta/r%Cbar 
         r%B = 1.0
         r%dM = -MW_g(i_CO2)                     ! 
         r%m = 1.0

         ! 2. Wet CO2 absorption        
         ! R2NH + CO2 + H2O_s <--> R2NH2 + HCO3
         ! This correspond to equations 3 (forward) and 4 (reverse)
         ! equations in mfix.dat

         r => reac(2)

         r%P_ref = P0 * unit%Pa                  ! dynes/cc        
         r%nv    = nv0 * unit%mol/unit%m3        ! mol-sites/m3
         r%Cbar  = r%nv                          ! mol/cc - 1 mol-R2NH3/mol-site

         r%E    = REAC_E(2)/Rgas_SI              ! K  
         r%dH   = REAC_DH(2)/Rgas_SI             ! K
         r%dS   = REAC_DS(2)/Rgas_SI             
         r%zeta = 10**REAC_LOGZETA(2)/unit%Pa   

         r%A = r%zeta/r%Cbar
         r%B = 1.0
         r%dM = -MW_g(i_CO2)  
         r%m = 1.0                    ! Chris Montgomery - initializzed exponent for loop below

         ! 3. Water physical absorption
         ! H2O(g) <=> H2O(phys)
         ! This correspond to equations 5 (forward) and 6 (reverse)
         ! equations in mfix.dat

         r => reac(3)

         r%P_ref = P0 * unit%Pa              ! dynes/cc        
         r%nv   = nv0 * unit%mol/unit%m3     ! mol-sites/m3
         r%Cbar = r%nv                       ! mol/cc - 1 mol-R2NH3/mol-site

         r%E    = REAC_E(3)/Rgas_SI              ! K  
         r%dH   = REAC_DH(3)/Rgas_SI             ! K
         r%dS   = REAC_DS(3)/Rgas_SI             
         r%zeta = 10**REAC_LOGZETA(3)/unit%Pa   

         r%A = r%zeta
         r%B = 1.0 * unit%mol/unit%m3
         r%dM = -MW_g(i_H2O_g)
         r%m = 1.0                    !Chris Montgomery - initializzed exponent for loop below

! Echo kinetic parameters on screen

         IF(myPE==PE_IO) THEN

            WRITE(*,*) " ===================================================== "
            WRITE(*,*) " Summary of kinetic parameters:"
            WRITE(*,*) ""
            WRITE(*,*) "    nv0 = ", nv0
            WRITE(*,*) ""
            WRITE(*,*) " 1. Dry CO2 absorption"
            WRITE(*,*) "    2R2NH + CO2 <--> R2NH2 + R2NCO2"
            WRITE(*,*) "    This correspond to equations 1 (forward) and 2 (reverse)"
            WRITE(*,*) "    equations in mfix.dat"
            WRITE(*,*) ""
            WRITE(*,*) "    E         = ",REAC_E(1)
            WRITE(*,*) "    DELTA H   = ",REAC_DH(1)
            WRITE(*,*) "    DELTA S   = ",REAC_DS(1)
            WRITE(*,*) "    LOG(ZETA) = ",REAC_LOGZETA(1)
            WRITE(*,*) ""
            WRITE(*,*) " 2. Wet CO2 absorption"
            WRITE(*,*) "    R2NH + CO2 + H2O_s <--> R2NH2 + HCO3"
            WRITE(*,*) "    This correspond to equations 3 (forward) and 4 (reverse)"
            WRITE(*,*) "    equations in mfix.dat"
            WRITE(*,*) ""
            WRITE(*,*) "    E         = ",REAC_E(2)
            WRITE(*,*) "    DELTA H   = ",REAC_DH(2)
            WRITE(*,*) "    DELTA S   = ",REAC_DS(2)
            WRITE(*,*) "    LOG(ZETA) = ",REAC_LOGZETA(2)
            WRITE(*,*) ""
            WRITE(*,*) " 3. Water physical absorption"
            WRITE(*,*) "    H2O(g) <=> H2O(phys)"
            WRITE(*,*) "    This correspond to equations 1 (forward) and 2 (reverse)"
            WRITE(*,*) "    equations in mfix.dat"
            WRITE(*,*) ""
            WRITE(*,*) "    E         = ",REAC_E(3)
            WRITE(*,*) "    DELTA H   = ",REAC_DH(3)
            WRITE(*,*) "    DELTA S   = ",REAC_DS(3)
            WRITE(*,*) "    LOG(ZETA) = ",REAC_LOGZETA(3)
            WRITE(*,*) " ===================================================== "

            IF(RUN_TYPE=='NEW') THEN ! SAVE VALUE OF nv0
               open (unit=888,file='nv0.dat')
               write(888,*)nv0
               close(888)
            ENDIF

         ENDIF


! Assigning Heats of Reactions to each (forward and backward) reaction   
! Although the user-defined HoR is not defined in mfix.dat,
! the calculation of HoR is turned off.
! The Solids phase receives the entire HoR
! Factor of 1000.used for unit conversion (kmole to mole)

      DO RR = 1, NO_OF_RXNS
         Reaction(RR)%calc_DH =  .FALSE.
         IF(.NOT.ALLOCATED(Reaction(RR)%HoR)) Allocate( Reaction(RR)%HoR(0:2))
         Reaction(RR)%HoR(0)  =  ZERO
         Reaction(RR)%HoR(2)  =  ZERO
      ENDDO

      Reaction(fwd_Dry_CO2_Adsorption)%HoR(1)  =  REAC_DH(1)*1000.0D0
      Reaction(rev_Dry_CO2_Adsorption)%HoR(1)  = -REAC_DH(1)*1000.0D0

      Reaction(fwd_Wet_CO2_Adsorption)%HoR(1)  =  REAC_DH(2)*1000.0D0
      Reaction(rev_Wet_CO2_Adsorption)%HoR(1)  = -REAC_DH(2)*1000.0D0

      Reaction(fwd_H2O_physisorption)%HoR(1)   =  REAC_DH(3)*1000.0D0
      Reaction(rev_H2O_physisorption)%HoR(1)   = -REAC_DH(3)*1000.0D0



    end subroutine 



    subroutine reaction_write()

    integer :: iReac

        do iReac = 1, nReac    
            write(*,*) iReac,reac(iReac)%A, reac(iReac)%E, reac(iReac)%dS, reac(iReac)%dH, reac(iReac)%dM
        end do

    end subroutine reaction_write








      END MODULE usr                                                                             
