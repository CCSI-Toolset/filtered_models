!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CONTINUUM_SOLIDS                                  !
!  Purpose: Check kinetic the run control namelist section             !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_CONTINUUM

! Global Variables:
!---------------------------------------------------------------------//
      USE constant
      USE run
      USE physprop

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: zero, undefined, undefined_i

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: LC
      DOUBLE PRECISION :: lsin_phi

!......................................................................!


! initialization of various dependent constants
      SIN_PHI = UNDEFINED   ! friction
      SIN2_PHI = UNDEFINED   ! schaeffer
      F_PHI = UNDEFINED    ! commented schaeffer section
      TAN_PHI_W = UNDEFINED   ! friction or jenkins

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_CONTINUUM")

! Check EP_star. This is used to populate ep_star_array which is what
! should be used elsewhere in the code. (see set_constprop)
      IF(EP_STAR == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'EP_STAR'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(EP_STAR < ZERO .OR. EP_STAR > ONE) THEN
         WRITE(ERR_MSG, 1001)'EP_STAR', iVal(EP_STAR)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! CHECK MU_s0
      IF (MU_S0 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MU_s0', MU_s0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! CHECK DIF_s0
      IF (DIF_S0 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'DIF_s0', DIF_s0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Constant solids viscosity is employed
      IF (MU_S0 .NE. UNDEFINED) THEN
! many of the solids phase transport coefficients that are needed for
! the granular energy equation will have be skipped in calc_mu_s. so
! do not evaluate granular energy if the solids viscosity is set to
! a constant.
         IF(GRANULAR_ENERGY) THEN
            WRITE(ERR_MSG,1100)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! needed by default solids-solids drag model
         IF (SMAX >=2) THEN
            IF (C_E == UNDEFINED) THEN
               WRITE(ERR_MSG,1101)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (C_F == UNDEFINED) THEN
               WRITE(ERR_MSG,1102)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ELSE
! PDE granular energy. Check kinetic theory specifications.
         IF (GRANULAR_ENERGY) THEN
            CALL CHECK_KT_TYPE
         ELSE
! Algebraic granular energy equation
            IF (C_E == UNDEFINED) THEN
               WRITE(ERR_MSG,1101)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
! needed by default solids-solids drag model. SMAX may be 1 for
! hybrid simulations and C_F is still needed.
            IF (SMAX >=2 .OR. DEM_SOLIDS) THEN
               IF (C_F == UNDEFINED) THEN
                  WRITE(ERR_MSG,1102)
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

 1100 FORMAT('Error 1100: Constant viscosity is specified but', /&
         'GRANULAR_ENERGY=.TRUE. Please correct the mfix.dat file')
 1101 FORMAT('Error 1101: Coefficient of restitution (C_E) not ',      &
         'specified.',/'Please correct the mfix.dat file.')
 1102 FORMAT('Error 1102: Coefficient of friction (C_F) not ',         &
         'specified.',/'Please correct the mfix.dat file.')
!

! If frictional stress modeling check various dependent/conflicting
! settings

! plastic/frictional stress model
      IF (FRICTION) THEN
         IF(SCHAEFFER) THEN
            WRITE(ERR_MSG, 1200)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
! Check that the granular energy PDE is solved.
         ELSEIF (.NOT.GRANULAR_ENERGY) THEN
            WRITE(ERR_MSG,1201)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
! Check the value specified for SAVAGE.
         ELSEIF(SAVAGE>2 .OR. SAVAGE<0) THEN
            WRITE(ERR_MSG, 1001)'SAVAGE', iVal(SAVAGE)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
! Verify that stress blending is not turned on.
         ELSEIF(BLENDING_STRESS) THEN
            WRITE(ERR_MSG, 1202)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(PHI == UNDEFINED) THEN
            WRITE(ERR_MSG, 1203)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
! used by friction bc
         ELSEIF(PHI_W == UNDEFINED) THEN
            WRITE(ERR_MSG, 1204)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! PHI & PHI_W are given in degrees but calculated in radian within
! the fortran codes
         SIN_PHI = SIN(PHI*PI/180.D0)
         TAN_PHI_W = TAN(PHI_W*PI/180.D0)
      ENDIF
 1201 FORMAT('Error 1201: The FRICTION solids stress model requires ', &
         /,'GRANULAR_ENERGY=.TRUE. Please correct the mfix.dat file.')
 1202 FORMAT('Error 1202: Cannot use BLENDING_STRESS with FRICTION ',&
         /,'Please correct the mfix.dat file.')
 1204 FORMAT('Error 1204: Angle of particle-wall friction (PHI_W) not',&
         ' specified.',/'Please correct the mfix.dat file.')


! plastic/frictional stress model
      IF(SCHAEFFER) THEN
         IF(FRICTION) THEN
            WRITE(ERR_MSG, 1200)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (PHI == UNDEFINED) THEN
            WRITE(ERR_MSG, 1203)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! PHI is given in degrees but calculated in radian within
! the fortran codes
         lsin_phi = sin(phi*PI/180.d0)
         SIN2_PHI = lSIN_PHI*lSIN_PHI
         F_PHI = (3.0D0 - 2.0D0*SIN2_PHI)/3.0D0    ! employed in commented
      ENDIF
 1200 FORMAT('Error 1200: FRICTION and SCHAEFFER models cannot be ',&
         'used',/'together. Please correct the mfix.dat file')
 1203 FORMAT('Error 1203: Angle of internal friction (PHI) not ',      &
         'specified.',/'Please correct the mfix.dat file.')


      IF(YU_STANDISH .AND. FEDORS_LANDEL) THEN
         WRITE(ERR_MSG, 1300)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(YU_STANDISH) THEN
! Yu_Standish correlation checks
         IF(SMAX < 2) THEN
            WRITE(ERR_MSG, 1301)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSEIF(FEDORS_LANDEL) THEN
! Fedors_Landel correlation checks.
         IF(SMAX /= 2) THEN
            WRITE(ERR_MSG, 1302)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF
 1300 FORMAT('Error 1300: FEDORS_LANDEL and YU_STANDISH correlations ',&
         'cannot be',/'used at the same time. Please correct the ',    &
         'mfix.dat file.')
 1301 FORMAT('Error 1301: YU_STANDISH correlation is for polydisperse',&
         ' mixtures',/'(MMAX >= 2). Please correct the mfix.dat file.')
 1302 FORMAT('Error 1302: FEDORS_LANDEL correlation is for binary ',   &
         'mixtures (MMAX=2).',/'Please correct the mfix.dat file.')


! Set the flags for blending stresses.
      IF(BLENDING_STRESS) THEN
! Turn off the default if SIGM_BLEND is set.
         IF(SIGM_BLEND)  TANH_BLEND = .FALSE.
      ELSE
         TANH_BLEND  = .FALSE.
         SIGM_BLEND  = .FALSE.
      ENDIF


      IF(MODEL_B) THEN
         DO LC = 1, MMAX
            IF(.NOT.CLOSE_PACKED(LC)) THEN
               WRITE(ERR_MSG, 1400) LC
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF
 1400 FORMAT('Error 1400: Solids phase ',I2,' is not CLOSE_PACKED.',/, &
         'All solids phases must be CLOSE_PACKED with MODEL_B=.TURE.',/ &
         'Please correct the mfix.dat file.')


! Check that phase number where added mass applies is properly defined.
      IF (ADDED_MASS) THEN
         IF(M_AM == UNDEFINED_I)THEN
            WRITE(ERR_MSG, 1500)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(M_AM == 0 .OR. M_AM > MMAX) THEN
            WRITE(ERR_MSG,1501)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1500 FORMAT('Error 1500: Must specify a disperse phase, M_AM, where ',&
         'the',/'virtual mass applies (ADDED_MASS).',/'Please correct',&
         ' the mfix.dat file.')
 1501 FORMAT('Error 1501: M_AM is out of range. [1,MMAX]',/'Please ',  &
         'correct the mfix.dat file.')
         ENDIF
      ENDIF


! Check name of radial distribution function
      SELECT CASE(trim(adjustl(RDF_TYPE)))

      CASE ('LEBOWITZ')
         RDF_TYPE_ENUM = LEBOWITZ

      CASE ('MODIFIED_LEBOWITZ')
         RDF_TYPE_ENUM = MODIFIED_LEBOWITZ

      CASE ('MANSOORI')
         RDF_TYPE_ENUM = MANSOORI

      CASE ('MODIFIED_MANSOORI')
         RDF_TYPE_ENUM = MODIFIED_MANSOORI

      CASE DEFAULT
            WRITE(ERR_MSG, 1600)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
       END SELECT

 1600 FORMAT('Error 1600: Unknown RDF_TYPE',/'Please ',  &
         'correct the mfix.dat file.')

! If the default (LEBOWITZ) is not set for a monodisperse case, then
! flag the error and exit. Otherwise, change it to CARNAHAN-STARLING.
      IF(MMAX == 1) THEN
         IF(RDF_TYPE_ENUM /= LEBOWITZ) THEN
            WRITE(ERR_MSG, 1601) trim(adjustl(RDF_TYPE))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1601 FORMAT('Error 1601: The RDF_TYPE should NOT be sepecified when ',&
         'MMAX = 1',/'because Carnahan-Starling is the only available',&
         ' radial distribution',/'function for monodisperse systems. ',&
         'Please correct the mfix.dat file.')

         ELSE
            RDF_TYPE_ENUM = CARNAHAN_STARLING
         ENDIF
      ENDIF

! Lane/Sarkar sub-grid models - 07/2015
      IF (SG_CYL_HYDRO) THEN

         IF ( MMAX > 2 ) THEN
            WRITE(ERR_MSG, 1700)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF ( SG_CYL_D == UNDEFINED ) THEN
            WRITE(ERR_MSG, 1000) 'SG_CYL_D'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE IF ( SG_CYL_D < ZERO ) THEN
            WRITE(ERR_MSG, 1001) 'SG_CYL_D', iVal(SG_CYL_D)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE IF ( UNITS == 'CGS' ) THEN
            SG_CYL_D = SG_CYL_D*0.1
         ENDIF

         IF ( SG_CYL_a == UNDEFINED ) THEN
            WRITE(ERR_MSG, 1000) 'SG_CYL_D'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE IF ( SG_CYL_a < ZERO ) THEN
            WRITE(ERR_MSG, 1001) 'SG_CYL_a', iVal(SG_CYL_a)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE IF ( UNITS == 'CGS' ) THEN
            SG_CYL_a= SG_CYL_a*0.1
         ENDIF

         IF ( SG_CYL_a <= 1.414*SG_CYL_D ) THEN
            WRITE(ERR_MSG, 1001) 'SG_CYL_a', iVal(SG_CYL_a)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         
         ! turn off the momentum and species flags for stationary cylinder porous phase
         MOMENTUM_X_EQ(2) = .FALSE.
         MOMENTUM_Y_EQ(2) = .FALSE.
         MOMENTUM_Z_EQ(2) = .FALSE.
         SPECIES_EQ(3) = .FALSE.

         !volume fraction occupied by the cylinders, should this be defined somewhere else??
         !constant during the simulation, throughout the domain, for now
         EP_cyl = 2*(PI*SG_CYL_D**2/4)/SG_CYL_a**2
         
      ENDIF 

      IF ( SG_CYL_ENERGY ) THEN
         IF ( .NOT.SG_CYL_HYDRO ) THEN
            WRITE(ERR_MSG, 1701)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF ( SG_CYL_T == UNDEFINED ) THEN
            WRITE(ERR_MSG, 1000) 'SG_CYL_T'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF ( SG_CYL_T < ZERO ) THEN
            WRITE(ERR_MSG, 1001) 'SG_CYL_T', iVal(SG_CYL_T)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF
      
 1700 FORMAT('Error 1700: MMAX must be equal to 2 for cylinder subgrid models.',/&
             'Please correct the mfix.dat file.')
 1701 FORMAT('Error 1701: SG_CYL_ENERGY requires SG_CYL_HYDRO = .TRUE.',/&
             'Please correct the mfix.dat file.')
! /Lane/Sarkar subgrid models

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')
 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG


      RETURN
      END SUBROUTINE CHECK_SOLIDS_CONTINUUM



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_KT_TYPE                                           !
!  Purpose: Check kinetic theory input specifications. These checks    !
!  are almost all related to the KT_TYPE keyword.                      !
!  Notes: To enter this routine granular_energy must be true           !
!                                                                      !
!  Author: J.Musser                                   Date: 04-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_KT_TYPE


! Global Variables:
!---------------------------------------------------------------------//
      USE constant
      USE run
      USE physprop

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: half, one, undefined

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! loop counters
      INTEGER :: I, J

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_KT_TYPE")

! These are some checks to satisfy legacy input:
      IF (AHMADI .AND. SIMONIN) THEN
         WRITE(ERR_MSG, 9001)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(AHMADI) THEN
         IF(KT_TYPE(1:6) /= 'AHMADI' .AND.                        &
            KT_TYPE(1:8) /= 'LUN_1984')THEN
            WRITE(ERR_MSG,9002)trim(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT = .TRUE.)
         ELSE
            KT_TYPE='AHMADI'
         ENDIF
      ELSEIF(SIMONIN) THEN
         IF(KT_TYPE(1:7) /= 'SIMONIN' .AND.                       &
            KT_TYPE(1:8) /= 'LUN_1984')THEN
            WRITE(ERR_MSG,9003)trim(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT = .TRUE.)
         ELSE
            KT_TYPE='SIMONIN'
         ENDIF
      ENDIF
 9001 FORMAT('Error 9001: Cannot specify AHMADI and SIMONIN together.',&
         /'Please correct the mfix.dat file.')
 9002 FORMAT('Error 9002: Cannot specify AHMADI and KT_TYPE = ',A,'.', &
         /'Please correct the mfix.dat file.')
 9003 FORMAT('Error 9003: Cannot specify SIMONIN and KT_TYPE = ',A,'.',&
         /'Please correct the mfix.dat file.')



! Check for valid options for kinetic theory models (KT_TYPE)
      SELECT CASE(trim(adjustl(KT_TYPE)))

!``````````````````````````````````````````````````````````````````````
      CASE ('IA_NONEP')
         KT_TYPE_ENUM = IA_2005
         IF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


!``````````````````````````````````````````````````````````````````````
      CASE ('GD_99')
         KT_TYPE_ENUM = GD_1999
         IF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(SMAX > 1) THEN
            WRITE(ERR_MSG,1002) TRIM(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


!``````````````````````````````````````````````````````````````````````
      CASE ('GTSH')
         KT_TYPE_ENUM = GTSH_2012
         IF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1002)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(SMAX > 1) THEN
            WRITE(ERR_MSG,1002) TRIM(KT_TYPE)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


!``````````````````````````````````````````````````````````````````````
      CASE ('GHD')
         KT_TYPE_ENUM = GHD_2007
! This variable is only used for GHD at this point...
! Define restitution coefficient matrix
         DO I = 1, SMAX
            DO J = 1, SMAX
               IF(r_p(I,J) == UNDEFINED) THEN
                  IF(C_E == UNDEFINED) THEN
                     WRITE(ERR_MSG,1003)
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ELSE
                     r_p(I,J) = C_e
                  ENDIF
               ENDIF
! just need to define r_p(1,2) and r_p(2,1) will be set.
               r_p(J,I) = r_p(I,J)
            ENDDO
         ENDDO

         IF(DRAG_TYPE_ENUM /= WEN_YU .AND. DRAG_TYPE_ENUM /= HYS) THEN
            WRITE(ERR_MSG, 1030)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(ADDED_MASS) THEN
            WRITE(ERR_MSG,1031)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(SMAX > 2) THEN  ! not sure this is still true!
            WRITE(ERR_MSG, 1032)
            CALL FLUSH_ERR_MSG
         ENDIF

! Automatically set SPECIES_EQ(MMAX) = .FALSE. to avoid potential
! checks/loops over the mmax species type eqn which has no meaning
         SPECIES_EQ(MMAX) = .FALSE.
         NMAX_s(MMAX) = 1

! currently set to avoid an overflow error in write_res0
! legacy variable?
         NMAX(MMAX) = 1

 1030 FORMAT('Error 1030: KT_TYPE = "GHD" is restricted to DRAG_TYPE', &
         'values of WEN_YU and HYS.',/'Please correct the mfix.dat ',  &
         'file.')
 1031 FORMAT('Error 1031: ADDED_MASS force cannot be applied with ',   &
         'GHD theory that',/'solves for mixture equations.',/'Please', &
         'correct the mifx.dat file.')
 1032 FORMAT('Warning 1032: GHD theory may not be valid for more ',    &
         'than two solids phases',/'it requires further development.')


!``````````````````````````````````````````````````````````````````````
      CASE ('AHMADI')
         KT_TYPE_ENUM = AHMADI_1995
         AHMADI = .TRUE.
         IF(.NOT.K_EPSILON) THEN
            WRITE(ERR_MSG,1040) 'K_EPSILON = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(C_F == UNDEFINED .AND. SMAX>=2) THEN
            WRITE(ERR_MSG, 1004)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1040 FORMAT('Error 1040: KT_TYPE = "AHMADI" requires ',A,/       &
         'Please correct the mfix.dat file.')


!``````````````````````````````````````````````````````````````````````
      CASE ('SIMONIN')
         KT_TYPE_ENUM = SIMONIN_1996
         SIMONIN = .TRUE.
         IF(.NOT.K_EPSILON) THEN
            WRITE(ERR_MSG,1050) 'K_EPSILON = .TRUE.'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(C_F == UNDEFINED .AND. SMAX>=2) THEN
            WRITE(ERR_MSG, 1004)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1050 FORMAT('Error 1050: KT_TYPE = "SIMONIN" requires ',A,/      &
         'Please correct the mfix.dat file.')


! Lun is the default implementation.
!``````````````````````````````````````````````````````````````````````
      CASE ('LUN_1984')
         KT_TYPE_ENUM = LUN_1984
! this version of the restitution coefficient is needed by most KT_TYPE
! models. it is also needed in default solids-solids drag model
         IF (C_E == UNDEFINED) THEN
            WRITE(ERR_MSG,1003)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(C_F == UNDEFINED .AND. SMAX>=2) THEN
            WRITE(ERR_MSG, 1004)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


      CASE DEFAULT
         WRITE(ERR_MSG,1001) trim(adjustl(KT_TYPE))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1001 FORMAT('Error 1001: Invalid or unknown KT_TYPE: ',A,/            &
         'Please correct the mfix.dat file.')

      END SELECT

! eventually this should be made specific to lun/ahmadi/simonin
! but because calc_gw_hw_cw in calc_u_friction is not consistent
! it currently must be defined for all kt_types whenever friction
! is invoked...
      ETA = (ONE + C_E)*HALF


 1002 FORMAT('Error 1002: KT_TYPE = ',A,' is for monodisperse',&
         ' solids',/'(MMAX = 1). Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Coefficient of restitution (C_E) not ',      &
         'specified.',/'Please correct the mfix.dat file.')

 1004 FORMAT('Error 1004: Coefficient of friction (C_F) not ',         &
         'specified.',/'Please correct the mfix.dat file.')


      RETURN
      END SUBROUTINE CHECK_KT_TYPE
