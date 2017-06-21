!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_gs                                                 C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Comments:                                                           C
!    If changes are made to this file, especially in regard to the     C
!    structure of any drag subroutine call, then such changes must     C
!    also be made to the analgous routine des_drag_gp in the file      C
!    des/drag_fgs.f for consistency!                                   C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DRAG_GS(M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar
      USE drag
      USE sendrecv
      USE discretelement
      USE ur_facs
      USE funits
      USE mms
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
! cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION :: USCM, VSCM, WSCM
! cell center value of x, y and z-particle velocity
      DOUBLE PRECISION :: USCM_HYS, VSCM_HYS, WSCM_HYS
! cell center value of U_g, V_g, W_g
      DOUBLE PRECISION :: UGC, VGC, WGC
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! drag coefficient
      DOUBLE PRECISION :: DgA
! current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
! indices of solids phases (continuous, discrete)
      INTEGER :: CM, DM, L
! temporary shift of total number of solids phases to account for both
! discrete and continuous solids phases used for the hybrid mdoel
      INTEGER :: MAXM
! tmp local variable for the particle diameter of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: DP_loc(2*DIM_M)
! tmp local variable for the solids volume fraction of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: EPs_loc(2*DIM_M)
! tmp local variable for the particle density of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: ROs_loc(2*DIM_M)
! correction factors for implementing polydisperse drag model
! proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION :: F_cor, tmp_sum, tmp_fac
! average particle diameter in polydisperse systems
      DOUBLE PRECISION :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION :: Y_i
! total solids volume fraction
      DOUBLE PRECISION :: phis
! temporary local variables to use for dummy arguments in subroutines
! void fraction, gas density, gas bulk density, solids volume fraction
! particle diameter, particle density
      DOUBLE PRECISION :: EPG, ROg, ROPg, EP_SM, DPM, ROs
!-----------------------------------------------

!!$omp  parallel do default(shared)                                   &
!!$omp  private( I,  IJK, IMJK, IJMK, IJKM, DM, MAXM, CM, L,          &
!!$omp           UGC, VGC, WGC, USCM, VSCM, WSCM, VREL, USCM_HYS,     &
!!$omp           VSCM_HYS, WSCM_HYS, tmp_sum, tmp_fac, Y_i, F_cor,    &
!!$omp           EP_SM, EPs_loc, ROs_loc, DP_loc, DPA, DPM, ROs,      &
!!$omp           phis, EPg, ROg, ROPg, Mu, DgA, F_gstmp)


      DO IJK = ijkstart3, ijkend3

         IF (FLUIDorP_FLOW_AT(IJK)) THEN

            I = I_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            MAXM = SMAX
            DO CM = 1,SMAX
               DP_loc(CM) = D_p(IJK,CM)
               EPs_loc(CM) = EP_S(IJK,CM)
               ROs_loc(CM) = RO_S(IJK,CM)
            ENDDO

! Append DEM phase information to the end for hybrid runs and update
! the number of solids 'phases' to loop over.
            IF(DES_CONTINUUM_HYBRID) THEN
               MAXM = SMAX + DES_MMAX
               DO DM = 1,DES_MMAX
                  L = SMAX + DM
                  DP_loc(L) = DES_D_p0(DM)
                  EPs_loc(L) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
                  ROs_loc(L) = DES_RO_S(DM)
               ENDDO
            ENDIF

! Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! magnitude of gas-solids relative velocity
            VREL = SQRT((UGC-USCM)**2 + (VGC-VSCM)**2 + (WGC-WSCM)**2)

! Laminar viscosity at a pressure boundary is given the value of the
! fluid cell next to it. This applies just to the calculation of the
! drag, in other routines the value of viscosity at a pressure boundary
! always has a zero value.
            IF (P_FLOW_AT(IJK)) THEN
               IF( FLUID_AT(EAST_OF(IJK)) ) THEN
                  Mu = MU_G(EAST_OF(IJK))
               ELSEIF ( FLUID_AT(WEST_OF(IJK)) ) THEN
                  Mu = MU_G(WEST_OF(IJK))
               ELSEIF ( FLUID_AT(NORTH_OF(IJK)) ) THEN
                  Mu = MU_G(NORTH_OF(IJK))
               ELSEIF ( FLUID_AT(SOUTH_OF(IJK)) ) THEN
                  Mu = MU_G(SOUTH_OF(IJK))
               ELSEIF ( FLUID_AT(TOP_OF(IJK)) ) THEN
                  Mu = MU_G(TOP_OF(IJK))
               ELSEIF ( FLUID_AT(BOTTOM_OF(IJK)) ) THEN
                  Mu = MU_G(BOTTOM_OF(IJK))
               ENDIF
            ELSE
               Mu = MU_G(IJK)
            ENDIF

! calculate the total solids volume fraction
            phis = ZERO
            DO L = 1, MAXM
! this is slightly /= one-ep_g due to round-off
               phis = phis + EPs_loc(L)
            ENDDO

! calculate the average paricle diameter and particle ratio
            DPA = ZERO
            tmp_sum = ZERO
            tmp_fac = ZERO
            DO L = 1, MAXM
               IF (phis .GT. ZERO) THEN
                  tmp_fac = EPs_loc(L)/phis
                  tmp_sum = tmp_sum + tmp_fac/DP_loc(L)
                ELSE
                  tmp_sum = tmp_sum + ONE/DP_loc(L) ! not important, but will avoid NaN's in empty cells
                ENDIF
            ENDDO
            DPA = ONE / tmp_sum
            Y_i = DP_loc(M)/DPA

! assign aliases for easy reference
            EPg = EP_G(IJK)
            ROg = RO_G(IJK)
            ROPg = ROP_G(IJK)
            EP_SM = EPs_loc(M)
            DPM = DP_loc(M)
            ROs = ROs_loc(M)


! determine the drag coefficient
            IF (EP_SM <= ZERO) THEN
               DgA = ZERO
            ELSEIF (EPg == ZERO) THEN
! this case will already be caught in most drag subroutines whenever
! RE==0 (for correlations in which RE includes EPg). however, this will
! prevent potential divisions by zero in some models by setting it now.
               DgA = ZERO
! Force a ZERO drag coefficient for MMS cases.
            ELSEIF (USE_MMS) THEN
               DgA = ZERO
            ELSE
               SELECT CASE(DRAG_TYPE_ENUM)

               CASE (SYAM_OBRIEN)
                  CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,DPM)

               CASE (GIDASPOW)
                  CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

               CASE (GIDASPOW_PCF)
                  CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)

               CASE (GIDASPOW_BLEND)
                  CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

               CASE (GIDASPOW_BLEND_PCF)
                  CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)

! Lane/Sarkar sub-grid models - 07/2015
               CASE (WEN_YU)
                  IF (SG_CYL_HYDRO) THEN
                     CALL DRAG_SG(DgA,EPg,MU,ROPg,VREL,DPM, &
                                  M,ROg,ROs,IJK)
                  ELSE
                     CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM,IJK)
                  ENDIF

               CASE (WEN_YU_PCF)
                  CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPA,IJK)

               CASE (KOCH_HILL)
                  CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)

               CASE (KOCH_HILL_PCF)
                  CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)

               CASE (BVK)
                  CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)

               CASE (USER_DRAG)
                  CALL DRAG_USR(IJK, M, DgA, EPg, Mu, ROg, VREL, DPM,  &
                     ROs, UGC, VGC, WGC)

               CASE (HYS)
! calculate velocity components of each solids phase
                  USCM_HYS = ZERO
                  VSCM_HYS = ZERO
                  WSCM_HYS = ZERO
                  IF(phis > ZERO) THEN
                     DO L = 1, MAXM
                        USCM_HYS = USCM_HYS + EPs_loc(L)*(UGC - &
                           AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I))
                        VSCM_HYS = VSCM_HYS + EPs_loc(L)*(VGC - &
                           AVG_Y_N(V_S(IJMK,L),V_S(IJK,L)))
                        WSCM_HYS = WSCM_HYS + EPs_loc(L)*(WGC - &
                           AVG_Z_T(W_S(IJKM,L),W_S(IJK,L)))
                     ENDDO
                     USCM_HYS = USCM_HYS/phis
                     VSCM_HYS = VSCM_HYS/phis
                     WSCM_HYS = WSCM_HYS/phis
                  ENDIF
! magnitude of gas-solids relative velocity
                  VREL = SQRT(USCM_HYS**2 +VSCM_HYS**2 +WSCM_HYS**2)

                  CALL DRAG_HYS(DgA,EPg,Mu,ROPg,VREL,&
                       DP_loc(:),DPA,Y_i,EPs_loc(:),phis,M,MAXM,IJK)

               CASE DEFAULT
                  CALL START_LOG
                  IF(.NOT.DMP_LOG) call open_pe_log(ier)
                  IF(DMP_LOG) WRITE (*, '(A,A)') &
                     'Unknown DRAG_TYPE: ', DRAG_TYPE
                  WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
                  CALL END_LOG
                  CALL mfix_exit(myPE)
               END SELECT   ! end selection of drag_type
            ENDIF   ! end if/elseif/else (ep_sm <= zero, ep_g==0)

! modify the drag coefficient to account for subgrid domain effects
            IF (SUBGRID_TYPE_ENUM .ne. UNDEFINED_SUBGRID_TYPE) THEN
               IF (SUBGRID_TYPE_ENUM .EQ. IGCI) THEN
                  CALL SUBGRID_DRAG_IGCI(DgA,EPg,Mu,ROg,DPM,ROs,IJK)
               ELSEIF (SUBGRID_TYPE_ENUM .EQ. MILIOLI) THEN
                  CALL SUBGRID_DRAG_MILIOLI(DgA,EPg,Mu,ROg,VREL,&
                       DPM,ROs,IJK)
               ENDIF
            ENDIF


! Modify drag coefficient to account for possible corrections and
! for differences between Model B and Model A
            IF(DRAG_TYPE_ENUM == HYS) THEN
! this drag model is handled differently than the others
               IF(Model_B)THEN
                  F_gstmp = DgA/EPg
               ELSE
                  F_gstmp = DgA
               ENDIF
            ELSE
               IF(DRAG_TYPE_ENUM == GIDASPOW_PCF .OR. &
                  DRAG_TYPE_ENUM == GIDASPOW_BLEND_PCF .OR. &
                  DRAG_TYPE_ENUM == WEN_YU_PCF .OR. &
                  DRAG_TYPE_ENUM == KOCH_HILL_PCF .OR. &
                  DRAG_TYPE_ENUM == BVK) THEN
! see erratum by Beetstra et al. (2007) : the correction factor differs
! for model A versus model B.
! application of the correction factor for model A is found from
! the correction factor for model B and neglects the Y_i**3 term
                  IF(Model_B) THEN
                     IF (M == 1) THEN
                        F_cor = (EPg*Y_i + phis*Y_i**2)
                     ELSE
                        F_cor = (EPg*Y_i + phis*Y_i**2 + &
                           0.064d0*EPg*Y_i**3)
                     ENDIF
                  ELSE
                     F_cor = Y_i
                  ENDIF
                  DgA = ONE/(Y_i*Y_i) * DgA * F_cor
               ENDIF

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g); all other models, eg., Wen_Yu
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_SM/EPg
               ELSE
                  F_gstmp = DgA * EP_SM                     !3D buoyancy model
               ENDIF
            ENDIF   !end if/else trim(drag_type=='hys')

! Determine drag force coefficient accounting for any under relaxation
            F_GS(IJK,M) = (ONE - UR_F_gs)*F_GS(IJK, M) + &
               UR_F_gs*F_gstmp

            IF(KT_TYPE_ENUM == GHD_2007) THEN
               IF(M==1) THEN
                  F_gs(IJK,MMAX) = F_gs(IJK,M)
               ELSE
                  F_gs(IJK,MMAX) = F_gs(IJK,MMAX) + F_gs(IJK,M)
               ENDIF
            ENDIF

         ELSE   ! .not.(fluidorp_flow_at(ijk)) branch

            F_GS(IJK,M) = ZERO
            IF(KT_TYPE_ENUM == GHD_2007) F_gs(IJK, MMAX) = ZERO

         ENDIF   ! end if (fluidorp_flow_at(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE DRAG_GS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function(s): C_DsxRe                                                C
!  Purpose:                                                            C
!     Calculate single sphere drag correlation multiplied by           C
!     the Reynolds number or                                           C
!     Calculate the single sphere drag correlation                     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

! Dalla Valle (1948)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DSXRE_DV(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

      C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
      RETURN
      END FUNCTION C_DSXRE_DV
!-----------------------------------------------------------------<<<

! Turton and Levenspiel (1986)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DSXRE_TL(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

      C_DSXRE_TL = 24.D0*(1.D0 + 0.173D0*RE**0.657D0) + &
         0.413D0*RE**2.09D0/(RE**1.09D0 + 16300.D0)
      RETURN
      END FUNCTION C_DSXRE_TL
!-----------------------------------------------------------------<<<


! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DS_SN(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

      C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
      RETURN
      END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SYAM_OBRIEN                                        C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Syamlal M, O'Brien TJ (1988). International Journal of           C
!        Multiphase Flow 14: 473-481.                                  C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SYAM_OBRIEN(lDgA,EPg,Mug,ROg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant, only : drag_c1, drag_d1, c
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: ldGA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1)
!     defined in the data file.
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variables which are function of EP_g
      DOUBLE PRECISION :: A, B
! Ratio of settling velocity of a multiparticle system to
! that of a single particle
      DOUBLE PRECISION :: V_rm
! Reynolds number
      DOUBLE PRECISION :: RE
!-----------------------------------------------
! External functions
!-----------------------------------------------
! Single sphere drag coefficient x Re
      DOUBLE PRECISION, EXTERNAL :: C_DSXRE_DV
!-----------------------------------------------

      IF(Mug > ZERO) THEN
         RE = DPM*VREL*ROg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! Calculate V_rm
      A = EPg**4.14D0
      IF (EPg <= 0.85D0) THEN
         B = drag_c1*EPg**1.28D0
      ELSE
        B = EPg**drag_d1
      ENDIF

      V_RM=HALF*(A-0.06D0*RE+&
           SQRT((3.6D-3)*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A) )

!------------------Begin cluster correction --------------------------
! uncomment the following four lines ...
!       V_RM = V_RM * (ONE + C(1)*&
!                      EXP(-A2*(RE-RE_C)**2 - A3*(EPg-EP_C)**2)* &
!                      RE*(1. - EPg))
!------------------End cluster correction ----------------------------

      lDgA = 0.75D0*Mug*EPg*C_DSXRE_DV(RE/V_RM)/(V_RM*DPM*DPM)

      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_SYAM_OBRIEN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW                                           C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.            C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW(lDgA,EPg,Mug,ROg,ROPg,VREL, DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM

!
      DOUBLE PRECISION, EXTERNAL :: C_DS_SN

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
!-----------------------------------------------


! Note the presence of gas volume fraction in ROPG
      RE = merge(DPM*VREL*ROPg/Mug, LARGE_NUMBER, MUg > ZERO)

! Dense phase
      IF(EPg <= 0.8D0) THEN
         lDgA = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
                1.75D0*ROg*VREL/DPM
      ELSE
! Dilute phase - EP_g >= 0.8
         IF(RE <= 1000D0)THEN
! this could be replaced with the function C_DS_SN
            C_d = C_DS_SN(RE)
         ELSE
            C_d = 0.44D0
         ENDIF
         lDgA = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM
      ENDIF

      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_GIDASPOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW_BLEND                                     C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!     original source unknown:                                         C
!     Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE  C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!     http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf. C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW_BLEND(lDgA,EPg,Mug,ROg,ROPg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant, only : PI
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
! Gidaspow switch function variables
      DOUBLE PRECISION :: Ergun
      DOUBLE PRECISION :: WenYu
      DOUBLE PRECISION :: PHI_gs
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! Dense phase - EP_g <= 0.8
      Ergun = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
               1.75D0*ROg*VREL/DPM
! Dilute phase - EP_g >= 0.8
      IF(RE <= 1000D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * &
           (ONE + 0.15D0 * RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF
      WenYu = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM

! Switch function
      PHI_gs = ATAN(150.D0*1.75D0*(EPg - 0.8D0))/PI + 0.5D0

! Blend the models
      lDgA = (1.D0-PHI_gs)*Ergun + PHI_gs*WenYu
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_GIDASPOW_BLEND



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_WEN_YU                                             C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium    C
!        Series 62: 100-111.                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_WEN_YU(lDgA,EPg,Mug,ROPg,VREL,DPM,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      ! subgrid models - Avik Sarkar - Jun 2014
      USE run
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE physprop
      USE rxns
      ! /subgrid models
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
! Subgrid models - Avik Sarkar June 2014
      INTEGER particlePhaseIdx, i, j, k, ijk, imjk, ijmk, ijkm
      DOUBLE PRECISION phis, h3DGasSolid_Sarkar
      DOUBLE PRECISION velTerminal, vSlipStar, filtSizeStar
      DOUBLE PRECISION hInf, alphaCoeff, VZero, powerVal
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      ! Sub-grid filtered model for 3D gas-particle flow - Avik Sarkar, July 2014
      h3DGasSolid_Sarkar = ZERO  ! by default, no subgrid correction
      
      IF(SG_3D_HYDRO) THEN

         i = I_OF(ijk)
         j = J_OF(ijk)
         k = K_OF(ijk)
         imjk = FUNIJK(i-1,j,k)
         ijmk = FUNIJK(i,j-1,k)
         ijkm = FUNIJK(i,j,k-1)         

         phis = ONE - EPg ! assuming only one solids phase
 
         velTerminal = ABS(GRAVITY)*DPM*DPM*(RO_S(IJK,1)-RO_g(IJK))/(18.0d0*MUg)
         vSlipStar = VREL/velTerminal

         filtSizeStar = DX(I_OF(ijk)) * DY(J_OF(ijk)) * DZ(K_OF(ijk))
         filtSizeStar = (2.0d0*filtSizeStar**(1.0d0/3.0d0))/(velTerminal**2/ABS(GRAVITY))

         hInf = 0.92 + 0.06497/(filtSizeStar**2)
         IF( hInf > 0.95 ) THEN
            hInf = 0.95
         ENDIF

         alphaCoeff = (0.7605 + 4.0560*phis - 13.6591*phis**2 + 16.1086*phis**3 - 3.1394*phis**4 )  &
              * ( ONE - exp(-300.0d0*phis) )  /  ( ONE + exp(100.0d0*(phis-0.55)) )  &
              * ( ONE - 0.5470/filtSizeStar + 0.2377/(filtSizeStar**2) )  &
              * ( ONE - 0.01811/(vSlipStar**2) )
         VZero = ( 0.3460 - 2.4062*phis ) / ( phis**0.1695 + 0.01 )  &
              * ( ONE + 0.2012/filtSizeStar - 0.08929/(filtSizeStar**2) )
         powerVal = ( 0.5380 - 0.2488*phis + 7.3176*phis**2 )  &
              * ( ONE + 0.4475/filtSizeStar - 0.4203/(filtSizeStar**2) )

         IF( vSlipStar < VZero ) THEN
            h3DGasSolid_Sarkar = ZERO 
         ELSE
            h3DGasSolid_Sarkar =  &
              hInf*(ONE - EXP( -alphaCoeff *(vSlipStar-VZero)**powerVal))
            IF(h3DGasSolid_Sarkar < ZERO) THEN 
               h3DGasSolid_Sarkar = ZERO 
            ENDIF
         ENDIF

      ENDIF
      ! /subgrid models

      IF(RE <= 1000.0D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * (ONE + 0.15D0*RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF

      lDgA = 0.75D0 * C_d * VREL * ROPg * EPg**(-2.65D0) / DPM
      IF (RE == ZERO) lDgA = ZERO

      ! subgrid model - AvS
      IF( SG_3D_HYDRO ) THEN
         lDgA = ( ONE - h3DGasSolid_Sarkar ) * lDgA
      ENDIF
      ! /subgrid models

      RETURN
      END SUBROUTINE DRAG_WEN_YU


! Lane/Sarkar sub-grid models - 07/2015
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SG                                                 C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SG(lDgA,EPg,MUg,ROPg,VREL,DPM,M,ROg,ROs,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! solids phase index
      INTEGER, INTENT(IN) :: M
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! solids density
      DOUBLE PRECISION, INTENT(IN) :: ROs
! cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
! Subgrid correction factor
      DOUBLE PRECISION :: SG_corr
! Ratio of solids
      DOUBLE PRECISION :: solidsRatio
! h2D correction
      DOUBLE PRECISION :: h2D
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      IF(RE <= 1000.0D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * (ONE + 0.15D0*RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF

      IF (EP_S(IJK,2) > ZERO) THEN  ! porous media present
! Sarkar subgrid models 2015
         SG_corr= ONE
         IF (M == 2) THEN
            SG_corr = ZERO
         ELSE
            solidsRatio = EP_S(IJK,1)/(ONE - EP_S(IJK,2)) !used repeatedly
            IF (solidsRatio > 0.0D0) then
               IF ((solidsRatio < 0.0012)) THEN
                  h2D = 2.7*(solidsRatio**0.234)
               ELSEIF ( solidsRatio < 0.014) THEN
                  h2D = -0.019*(solidsRatio**(-0.455)) + 0.963
               ELSEIF (solidsRatio < 0.25) THEN 
                  h2D = 0.868*exp(-0.38*solidsRatio) - 0.176*exp(-119.2*solidsRatio)
               ELSEIF (solidsRatio < 0.30) THEN 
                  h2D = -4.59e-5*exp(19.75*solidsRatio) + 0.852*exp(-0.268*solidsRatio)
               ELSEIF (solidsRatio <= 0.64) THEN
                  h2D = (-0.4341*solidsRatio + 0.8998) * (1 - exp(42.68*(solidsRatio-0.64)))
               ELSE 
                  h2D = ZERO
               ENDIF
            ELSE
               h2D = ZERO
            ENDIF
            SG_corr = (ONE - EP_S(IJK,1))*(ONE - h2D)
         ENDIF

      ELSE
         CALL SUBGRID_DRAG_IGCI(lDgA,EPg,Mug,ROg,DPM,ROs,IJK)
         RETURN
      ENDIF

      lDgA = 0.75D0 * C_d * VREL * ROPg * EPg**(-2.65D0) / DPM
      IF (RE == ZERO) lDgA = ZERO

      ! If sub-grid filtered model is used, multiply by the correction factor
      IF (EP_S(IJK,2) > ZERO) THEN
         lDgA = SG_corr*lDgA
      ENDIF

      RETURN
      END SUBROUTINE DRAG_SG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_IGCI                                       C
!  Purpose: Calculate subgrid correction to the gas-solids drag        C
!           coefficient developed by Wen-Yu                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes, e.g., fix inconsistenty with analogous      C
!     calls in des_drag_gp and with new variable density feature       C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., Pannala, S., Benyahia, S., & Sundaresan S.,            C
!        Validation studies on filtered model equations for gas-       C
!        particle flows in risers, Industrial & Engineering Chemistry  C
!        Research, 2012, 51(4), 2094-2103                              C
!                                                                      C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_IGCI(lDgA,EPg,Mug,ROg,DPM,ROs,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run, only : filter_size_ratio, SUBGRID_WALL
      USE constant, only : GRAVITY
      USE geometry, only : VOL,AXY,DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(INOUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects
      DOUBLE PRECISION :: F_Subgrid
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION :: F_SubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION :: vt
! filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
! inverse Froude number, or dimensionless filter size
      DOUBLE PRECISION :: Inv_Froude
! total solids volume fraction
      DOUBLE PRECISION :: EPs
! Variables for Igci model
      DOUBLE PRECISION :: GG_phip, h_phip, h_phip2, c_function,&
                          f_filter
!-----------------------------------------------

! initialize
      F_Subgrid = ONE
      F_SubgridWall = ONE

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
      vt = GRAVITY*DPM*DPM*(ROs - ROg) / (18.0d0*Mug)
! filter size calculation for each specific gridcell volume
      IF(DO_K) THEN
         filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
      ELSE
         filtersize = filter_size_ratio * DSQRT(AXY(IJK))
      ENDIF

! dimensionless inverse of Froude number
      IF(ABS(vt) > SMALL_NUMBER) THEN
         Inv_Froude =  filtersize * GRAVITY / vt**2
      ELSE
         Inv_Froude =  LARGE_NUMBER
      ENDIF

! total solids volume fraction
      EPs = ONE - EPg

      IF (EPs .LT. 0.0012d0) THEN
         h_phip = 2.7d0*(EPs**0.234)
      ELSEIF (EPs .LT. 0.014d0) THEN
         h_phip = -0.019d0/(EPs**0.455) + 0.963d0
      ELSEIF (EPs .LT. 0.25d0) THEN
         h_phip = 0.868d0*EXP((-0.38*EPs)) - &
            0.176d0*EXP((-119.2*EPs))
      ELSEIF (EPs .LT. 0.455d0) THEN
         h_phip = -4.59d-5*EXP((19.75*EPs)) + &
            0.852d0*EXP((-0.268*EPs))
      ELSEIF (EPs .LE. 0.59d0) THEN
         h_phip = (EPs - 0.59d0) * (-1501.d0*(EPs**3) + &
            2203.d0*(EPs**2) - 1054.d0*EPs + 162.d0)
      ELSE
         h_phip=ZERO
      ENDIF

      IF (EPs .LT. 0.18d0) THEN
          GG_phip = (EPs**0.24)*(1.48d0 + EXP(-18.0*EPs))
      ELSE
          GG_phip = ONE
      ENDIF

! a filter function needed in Igci Filtered/subgrid Model [dimensionless]
      f_filter = (Inv_Froude**1.6) / ((Inv_Froude**1.6)+0.4d0)
      h_phip2=h_phip*GG_phip
      c_function=-h_phip2*f_filter
      F_Subgrid =(ONE + c_function)

      IF (SUBGRID_WALL) THEN
         CALL SUBGRID_DRAG_WALL(F_SubgridWall,vt,IJK)
      ENDIF

      lDgA = F_SubgridWall*F_Subgrid * lDgA

      RETURN
      END SUBROUTINE SUBGRID_DRAG_IGCI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_MILIOLI                                    C
!  Purpose: Calculate subgrid correction to the gas-solids drag        C
!           coefficient developed by Wen-Yu                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes, e.g., fix inconsistenty with analogous      C
!     calls in des_drag_gp and with new variable density feature       C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Milioli, C. C., et al., Filtered two-fluid models of fluidized   C
!        gas-particle flows: new constitutive relations, AICHE J,      C
!        doi: 10.1002/aic.14130                                        C
!                                                                      C
!  Comments:                                                           C
!     Still needs to be reviewed for accuracy with source material     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_MILIOLI(lDgA,EPg,Mug,ROg,VREL,DPM,ROs,&
         IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run, only : filter_size_ratio, SUBGRID_WALL
      USE constant, only : GRAVITY
      USE geometry, only : VOL,AXY,DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(INOUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects
      DOUBLE PRECISION :: F_Subgrid
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION :: F_SubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION :: vt
! filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
! inverse Froude number, or dimensionless filter size
      DOUBLE PRECISION :: Inv_Froude
! dimensionless slip velocity = VREL/vt
      DOUBLE PRECISION :: vslip
! total solids volume fraction
      DOUBLE PRECISION :: EPs
! Variables for Milioli model
      DOUBLE PRECISION :: h1, henv, hlin
!-----------------------------------------------

! initialize
      F_Subgrid = ONE
      F_SubgridWall = ONE

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
      vt = GRAVITY*DPM*DPM*(ROs - ROg) / (18.0d0*Mug)
! filter size calculation for each specific gridcell volume
      IF(DO_K) THEN
         filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
      ELSE
         filtersize = filter_size_ratio * DSQRT(AXY(IJK))
      ENDIF
! dimensionless inverse of Froude number
      IF(ABS(vt) > SMALL_NUMBER) THEN
         Inv_Froude =  filtersize * GRAVITY / vt**2
      ELSE
         Inv_Froude =  LARGE_NUMBER
      ENDIF
! total solids volume fractionn
      EPs = ONE - EPg
! dimensionless slip velocity between gas and solids phase M
      Vslip = VREL / vt

      IF (Inv_Froude .LE. 1.028d0) THEN
         h1 = (1.076d0 + 0.12d0*Vslip - (0.02d0/(Vslip+0.01d0)))*EPs + &
            (0.084d0 + 0.09d0*Vslip - (0.01d0/(0.1d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.53d0) THEN
            henv = (6.8d0*(ONE+EPs)*(EPs**0.3)) / &
               (10.d0*(EPs**1.5) + 5.d0)
         ELSEIF (EPs .GT. 0.53d0 .AND. EPs .LE. 0.65d0) THEN
            henv = (2.23d0*((0.65d0-EPs)**(0.45))) / &
               ((ONE/EPs)-ONE)
         ELSEIF (EPs .GT. 0.65d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 1.028d0 .AND. &
              Inv_Froude .LE. 2.056d0) THEN
         h1 = (1.268d0 - (0.2d0*Vslip) + (0.14d0/(Vslip+0.01d0)))*EPs + &
            (0.385d0 + 0.09d0*Vslip - (0.05d0/(0.2d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.53d0) THEN
            henv = (8.6d0*(ONE+EPs)*(EPs**0.2)) / (10.d0*EPs + 6.3d0)
         ELSEIF (EPs .GT. 0.53d0 .AND. EPs .LE. 0.65d0) THEN
            henv = (0.423d0*((0.65d0-EPs)**0.3)) / (ONE-(EPs**0.4))
         ELSEIF (EPs .GT. 0.65d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 2.056d0 .AND. &
              Inv_Froude .LE. 4.112d0) THEN
         h1 = ((0.018d0*Vslip + 0.1d0)/(0.14d0*Vslip + 0.01d0))*EPs + &
            (0.9454d0 - (0.09d0/(0.2d0*Vslip + 0.01d0)))
         IF (EPs .LE. 0.5d0) THEN
            henv = (7.9d0*(ONE+EPs)*(EPs**0.2)) / &
               (10.d0*(EPs**0.9) + 5.d0)
         ELSEIF (EPs .GT. 0.5d0 .AND. EPs .LE. 0.63d0) THEN
            henv = (0.705d0*((0.63d0-EPs)**0.3)) / (ONE-(EPs**0.7))
         ELSEIF (EPs .GT. 0.63d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 4.112d0 .AND. &
              Inv_Froude .LE. 8.224d0) THEN
         h1 = ((0.05d0*Vslip+0.3d0)/(0.4d0*Vslip+0.06d0))*EPs + &
            (0.9466d0 - (0.05d0/(0.11d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.45d0) THEN
            henv = (7.9d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.6)) + 3.6d0)
         ELSEIF (EPs .GT. 0.45d0 .AND. EPs .LE. 0.57d0) THEN
            henv = (0.78d0*((0.57d0-EPs)**0.2)) / (ONE-(EPs**0.9))
         ELSEIF (EPs .GT. 0.57d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 8.224d0 .AND. &
              Inv_Froude .LE. 12.336d0) THEN
         h1 = ((1.3d0*Vslip+2.2d0)/(5.2d0*Vslip+0.07d0))*EPs + &
            (0.9363d0-(0.11d0/(0.3d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.35d0) THEN
            henv = (7.6d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.6)) + 3.3d0)
         ELSEIF (EPs .GT. 0.35d0 .AND. EPs .LE. 0.55d0) THEN
            henv = (0.81d0*((0.55d0-EPs)**0.3)) / (ONE-(EPs**0.7))
         ELSEIF (EPs .GT. 0.55d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 12.336d0 .AND. &
              Inv_Froude .LE. 16.448d0) THEN
         h1 = ((2.6d0*Vslip+4.d0)/(10.d0*Vslip+0.08d0))*EPs + &
            (0.926d0-(0.17d0/(0.5d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (8.4d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.5)) + 3.3d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (1.01d0*((0.52d0-EPs)**0.03))/(ONE-(EPs**0.9))
         ELSEIF (EPs .GT. 0.52d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 16.448d0 .AND. &
              Inv_Froude .LE. 20.56d0) THEN
         h1 = ((2.5d0*Vslip+4.d0)/(10.d0*Vslip+0.08d0))*EPs + &
            (0.9261d0-(0.17d0/(0.5d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (8.4d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.5)) + 3.3d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (1.065d0*((0.52d0-EPs)**0.3))/(ONE-EPs)
         ELSEIF (EPs .GT.0.52d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 20.56d0) THEN
         h1 = ((1.6d0*Vslip+4.d0)/(7.9d0*Vslip+0.08d0))*EPs + &
            (0.9394d0 - (0.22d0/(0.6d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (9.d0*(ONE+EPs)*(EPs**0.15)) / &
               (10.d0*(EPs**0.45) + 4.2d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (0.91d0*((0.52d0-EPs)**0.4))/(ONE-(EPs**0.6))
         ELSEIF (EPs .GT. 0.52d0) THEN
            henv=ZERO
         ENDIF
      ENDIF

      IF (h1 .GT. ZERO) THEN
         hlin=h1
      ELSE
         hlin=ZERO
      ENDIF

      IF (Inv_Froude .LT. 1.028d0) THEN
! for very small filtered size, the drag wont be changed:
! F_Subgrid = 1.0 - H where H = 0.0
         F_Subgrid = ONE
      ELSE
! MIN(henv,hlin) is H in Milioli paper, 2013
         F_Subgrid = ONE - MIN(henv,hlin)
      ENDIF

! Filtered drag = (1 - H)*Microscopic_drag; it is strange Milioli takes EPs
!     F_Subgrid = EPs*(ONE-hmili)

      IF (SUBGRID_WALL) THEN
         CALL SUBGRID_DRAG_WALL(F_SubgridWall,vt,IJK)
      ENDIF

      lDgA = F_SubgridWall*F_Subgrid * lDgA


      RETURN
      END SUBROUTINE SUBGRID_DRAG_MILIOLI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_WALL                                       C
!  Purpose: Calculate subgrid correction arising from wall to the      C
!     gas-solids drag coefficient                                      C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., and Sundaresan, S., Verification of filtered two-      C
!        fluid models for gas-particle flows in risers, AICHE J.,      C
!        2011, 57 (10), 2691-2707.                                     C
!                                                                      C
!  Comments: Currently only valid for free-slip wall but no checks     C
!     are made to ensure user has selected free-slip wall when this    C
!     option is invoked                                                C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_WALL(lSubgridWall,vt,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant, only : GRAVITY
      USE cutcell, only : DWALL
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION, INTENT(OUT) :: lSubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION, INTENT(IN) :: vt
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! values are only correct for FREE-Slip walls
      DOUBLE PRECISION, PARAMETER :: a22=6.0d0, b22=0.295d0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! dimensionless distance to wall
      DOUBLE PRECISION :: x_d
!-----------------------------------------------
! initialize
      lSubgridWall = ONE

! dimensionless distance to the Wall
      x_d = DWALL(IJK) * GRAVITY / vt**2

! decrease exponentionally away from the wall
! more complex model could be implemented with JJ wall model
      lSubgridWall = ONE / ( ONE + a22 * (EXP(-b22*x_d)) )

      RETURN
      END SUBROUTINE SUBGRID_DRAG_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_KOCH_HILL                                          C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Clay Sutton (Lehigh University)            Date: 14-Jul-04  C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Sofiane Benyahia                           Date: 21-Jan-05  C
!                                                                      C
!  Literature/Document References:                                     C
!     Benyahia S, Syamlal M, O'Brien TJ (2006). Powder Technology      C
!        162: 166-174.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 213-241.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 243-278.                                                 C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_KOCH_HILL(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter if pcf otherwise DPM again
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases
      DOUBLE PRECISION, INTENT(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! transition Reynolds numbers
      DOUBLE PRECISION :: Re_Trans_1, Re_Trans_2
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! zero Re function for low Reynolds number
      DOUBLE PRECISION :: F_0
! inertial function for low Reynolds number
      DOUBLE PRECISION :: F_1
! zero Re function for high Reynolds number
      DOUBLE PRECISION :: F_2
! inertial function for high Reynolds number
      DOUBLE PRECISION :: F_3
! dimensionless drag force F
      DOUBLE PRECISION :: F
! weighting factor to compute F_0 and F_2
      DOUBLE PRECISION :: w
!-----------------------------------------------


      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG and factor of 1/2
         RE = 0.5D0*DPA*VREL*ROPG/Mug        ! if pcf DPA otherwise DPM
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      F_STOKES = 18.D0*Mug*EPg*EPg/DPM**2    ! use DPM
      w = EXP(-10.0D0*(0.4D0-phis)/phis)

      IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
         F_0 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
            135.0D0/64.0D0*phis*LOG(phis) + 17.14D0*phis) / &
            (1.0D0 + 0.681D0*phis - 8.48D0*phis*phis + &
            8.16D0*phis**3) + w*10.0D0*phis/(1.0D0-phis)**3
      ELSEIF(phis >= 0.4D0) THEN
         F_0 = 10.0D0*phis/(1.0D0-phis)**3
      ENDIF

      IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
        F_1 = dsqrt(2.0D0/phis) / 40.0D0
      ELSE IF(phis > 0.1D0) THEN
        F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
      ENDIF

      IF(phis < 0.4D0) THEN
        F_2 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
           135.0D0/64.0D0*phis*LOG(phis) + 17.89D0*phis) / &
           (1.0D0 + 0.681D0*phis - 11.03D0*phis*phis + &
           15.41D0*phis**3)+ w*10.0D0*phis/(1.0D0-phis)**3
      ELSE
         F_2 = 10.0D0*phis/(1.0D0-phis)**3
      ENDIF

      IF(phis < 0.0953D0) THEN
         F_3 = 0.9351D0*phis + 0.03667D0
      ELSE
         F_3 = 0.0673D0 + 0.212D0*phis +0.0232D0/(1.0-phis)**5
      ENDIF

      Re_Trans_1 = (F_2 - 1.0D0)/(3.0D0/8.0D0 - F_3)
      Re_Trans_2 = (F_3 + dsqrt(F_3*F_3 - 4.0D0*F_1 &
           *(F_0-F_2))) / (2.0D0*F_1)

      IF(phis <= 0.01D0 .AND. RE <= Re_Trans_1) THEN
         F = 1.0D0 + 3.0D0/8.0D0*RE
      ELSEIF(phis > 0.01D0 .AND. RE <= Re_Trans_2) THEN
         F = F_0 + F_1*RE*RE
      ELSEIF(phis <= 0.01D0 .AND. RE > Re_Trans_1 .OR.   &
         phis >  0.01D0 .AND. RE > Re_Trans_2) THEN
         F = F_2 + F_3*RE
      ELSE
         F = zero
      ENDIF

      lDgA = F * F_STOKES
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_KOCH_HILL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_BVK                                                C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62           C
!     (Jan 2007)                                                       C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_BVK(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases
      DOUBLE PRECISION, INTENT(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! dimensionless drag force F
      DOUBLE PRECISION :: F
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPA*VREL*ROPg/Mug        ! use DPA
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! eq(9) BVK J. fluid. Mech. 528, 2005
! (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
      F_STOKES = 18D0*Mug*EPg/DPM**2   ! use DPM

      F = 10d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))
      F = F + 0.413d0*RE/(24.d0*EPg**2) * &
             (ONE/EPg + 3d0*EPg*phis + 8.4d0/RE**0.343) / &
             (ONE+10.d0**(3d0*phis)/RE**(0.5+2.d0*phis))

      lDgA = F*F_STOKES
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_BVK



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_HYS                                                C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Yin, X, Sundaresan, S. (2009). AIChE Journal 55: no 6, 1352-     C
!     1368                                                             C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_HYS(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,Y_I,EP_sM,PHIS,M,MAXM,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE drag, only : beta_ij
      USE run, only : LAM_HYS
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! local variable for the particle diameter
      DOUBLE PRECISION :: DPM(2*DIM_M)
! average particle diameter
      DOUBLE PRECISION, INTENT(IN) :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION, INTENT(IN) :: Y_i
! local variable for the solids volume fraction
      DOUBLE PRECISION :: EP_sM(2*DIM_M)
! total solids volume fraction of solids phases
      DOUBLE PRECISION, INTENT(IN) :: PHIS
! current solids phase index and fluid cell index
      INTEGER, INTENT(IN) :: M, IJK
! maximum number of solids phases
      INTEGER, INTENT(IN) :: MAXM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! minimum particle diameter in mixture
      DOUBLE PRECISION :: DPmin
! Index for particles of other species
      INTEGER :: L
! Reynolds number
      DOUBLE PRECISION :: RE
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! dimensionless drag force F
      DOUBLE PRECISION :: F
! Polydisperse correction factor for YS drag relation
      DOUBLE PRECISION :: a_YS
! Lubrication interaction prefactor in YS drag relation
      DOUBLE PRECISION :: alpha_YS
! Friction coefficient for a particle of type i (HYS drag relation)
      DOUBLE PRECISION :: beta_i_HYS
! Friction coefficient for a particle of type j (HYS drag relation)
      DOUBLE PRECISION :: beta_j_HYS
! Stokes drag of a particle of type j
      DOUBLE PRECISION :: FSTOKES_j
! Diameter ratio for particle of type j
      DOUBLE PRECISION :: Y_i_J
! Variable for Beetstra et. al. drag relation
      DOUBLE PRECISION :: F_D_BVK
! Variable for YS drag relation
      DOUBLE PRECISION :: F_YS
!-----------------------------------------------

      IF (Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPA*VREL*ROPg/Mug   ! use DPA
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! (this F_Stokes is /= of Koch_Hill by a factor of ep_g/ep_sm)
! (this F_Stokes is /= of BVK by a factor of 1/ep_sm)
      F_STOKES = 18D0*Mug*EPg*EP_SM(M)/DPM(M)**2   ! use DPM

! Find smallest diameter if number of particle types is greater than 1
      Dpmin= DPM(1)
      IF (MAXM > 1) THEN
         DO L=2,MAXM
            Dpmin = MIN(Dpmin,DPM(L))
         ENDDO
      ENDIF

      a_YS = 1d0 - 2.66d0*phis + 9.096d0*phis**2 - 11.338d0*phis**3

! Calculate the prefactor of the off-diagonal friction coefficient
! Use default value of lamdba if there are no particle asparities
      alpha_YS = 1.313d0*LOG10(DPmin/lam_HYS) - 1.249d0


! Beetstra correction for monodisperse drag
      F_D_BVK = ZERO
      F = 10.d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))

      IF(RE > ZERO) F_D_BVK = F + 0.413d0*RE/(24.d0*EPg**2)*&
         (ONE/EPg + 3.d0*EPg*phis + 8.4d0/RE**0.343d0) / &
         (ONE+10**(3.d0*phis)/RE**(0.5d0+2.d0*phis))

! YS correction for polydisperse drag
      F_YS = 1d0/EPg + (F_D_BVK - 1d0/EPg)*&
                       (a_YS*Y_i+(1d0-a_YS)*Y_i**2)
      F_YS = F_YS*F_STOKES
      beta_i_HYS = F_YS

      DO L= 1,MAXM
         IF (L /= M) THEN
            Y_i_J = DPM(L)/DPA
            beta_j_HYS = 1.d0/EPg + (F_D_BVK - 1.d0/EPg) * &
               (a_YS*Y_i_J + (1d0-a_YS)*Y_i_J**2)
            FSTOKES_j = 18.D0*Mug*EP_sM(L)*EPg/&
               DPM(L)**2

            beta_j_HYS = beta_j_HYS*FSTOKES_j

! Calculate off-diagonal friction coefficient
            beta_ij(IJK,M,L) = ZERO

! This if statement prevents NaN values from appearing for beta_ij
            IF (EP_sM(M) > ZERO .AND. EP_SM(L) > ZERO) &
               beta_ij(IJK,M,L) = (2.d0*alpha_YS*EP_sM(M)*EP_sM(L))/ &
                  (EP_sM(M)/beta_i_HYS + EP_sM(L)/beta_j_HYS)
            F_YS = F_YS + beta_ij(IJK,M,L)

         ENDIF   ! end if (J/=M)
      ENDDO   ! end do (J=1,MAXM)


      lDgA = F_YS

      RETURN
      END SUBROUTINE DRAG_HYS
