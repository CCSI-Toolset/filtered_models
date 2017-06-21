!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_ENERGY_EQ                                         C
!  Purpose: Solve energy equations                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Eliminate energy calculations when doing DEM               C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_ENERGY_EQ(IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs
      USE pgcor
      USE pscor
      USE leqsol
      USE bc
      USE energy
      USE rxns
      Use ambm
      Use tmp_array, S_p => ARRAY1, S_C => ARRAY2, EPs => ARRAY3
      Use tmp_array1, VxGama => ARRAYm1
      USE compar
      USE discretelement
      USE des_thermo
      USE mflux
      USE mpi_utility
      USE sendrecv
      USE ps
      USE mms
      USE functions
      USE constant  ! Lane/Sarkar sub-grid models - 07/2015

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M
      INTEGER :: TMP_SMAX
!  Cp * Flux
      DOUBLE PRECISION :: CpxFlux_E(DIMENSION_3), &
                          CpxFlux_N(DIMENSION_3), &
                          CpxFlux_T(DIMENSION_3)
! previous time step term
      DOUBLE PRECISION :: apo
! Indices
      INTEGER :: IJK, I, J, K
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! Arrays for storing errors:
! 120 - Gas phase energy equation diverged
! 121 - Solids energy equation diverged
! 12x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! Lane sub-grid models - 07/2015
      DOUBLE PRECISION :: EP_G_STAR, V_G_MAG, RE_CG, PR_CG, NU_CG, H_CG
      DOUBLE PRECISION :: L_R, A_C
      DOUBLE PRECISION :: Q_CG

      DOUBLE PRECISION :: L_STAR
      DOUBLE PRECISION :: EP_S_STAR, V_S_STAR, PE_CS, NU_CS, H_CS
      DOUBLE PRECISION :: T_SUSP, Q_CS

! temporary use of global arrays:
! arraym1 (locally vxgama)
! the volume x average gas-solids heat transfer at cell centers
!      DOUBLE PRECISION :: VXGAMA(DIMENSION_3, DIMENSION_M)
! array1 (locally s_p)
! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
!      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c)
! source vector: constant part becomes part of b_m vector
!      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
!      DOUBLE PRECISION :: EPS(DIMENSION_3)

! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      call lock_ambm         ! locks arrys a_m and b_m
      call lock_tmp_array    ! locks arraym1 (locally vxgama)
      call lock_tmp_array1   ! locks array1, array2, array3
                             ! (locally s_p, s_c, eps)
! Initialize error flags.
      Err_l = 0

      TMP_SMAX = SMAX
      IF(DISCRETE_ELEMENT) THEN
         TMP_SMAX = 0   ! Only the gas calculations are needed
      ENDIF

! Lane/Sarkar sub-grid models - 08/2015
      IF (SG_CYL_ENERGY) THEN
         TMP_SMAX = 1
      ENDIF

! initializing
      DO M = 0, TMP_SMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER)
      ENDDO

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gE(IJK)
            ELSE
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gSE(IJK)
            ENDIF
         ENDIF

         IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gN(IJK)
            ELSE
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gSN(IJK)
            ENDIF
         ENDIF

         IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gT(IJK)
            ELSE
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gST(IJK)
            ENDIF
         ENDIF

         IF (FLUID_AT(IJK)) THEN
            APO = ROP_GO(IJK)*C_PG(IJK)*VOL(IJK)*ODT
            S_P(IJK) = APO + S_RPG(IJK)*VOL(IJK)
            S_C(IJK) = APO*T_GO(IJK)-HOR_G(IJK)*VOL(IJK)+S_RCG(IJK)*VOL(IJK)
! Lane sub-grid models - 07/2015
            IF (SG_CYL_ENERGY) THEN
               IF (EP_S(IJK,1) <= DIL_EP_S .AND. EP_S(IJK,2) > ZERO) THEN
                  ! Variables for Nusselt number calculation
                  EP_G_STAR = EP_G(IJK)/(1.0-EP_S(IJK,2))
                  V_G_MAG = SQRT(U_G(IJK)**2 + V_G(IJK)**2 + W_G(IJK)**2)
                  RE_CG = RO_g(IJK)*V_G_MAG*SG_CYL_D/MU_G(IJK)
                  PR_CG = C_PG(IJK)*MU_G(IJK)/K_G(IJK)
                  ! Nusselt number and heat transfer coefficient
                  NU_CG = (0.289 + 2.53*(SG_CYL_D/SG_CYL_A)**1.65) * RE_CG**0.564 * PR_CG**0.430
                  H_CG = NU_CG*K_G(IJK)/SG_CYL_D
                  ! Area calculation and correction
                  L_R = DX(I)*DY(J)/SG_CYL_A**2  ! Length-scaling ratio
                  A_C = 2.0*PI*(L_R*SG_CYL_D)  ! 2D/3D Correction factor
                  ! Source term
                  Q_CG = EP_G_STAR*H_CG*A_C*(SG_CYL_T-T_G(IJK))
                  S_C(IJK) = S_C(IJK) + Q_CG
               ENDIF
            ENDIF
            IF(USE_MMS) S_C(IJK) = S_C(IJK) + MMS_T_G_SRC(IJK)*VOL(IJK)
         ELSE
            S_P(IJK) = ZERO
            S_C(IJK) = ZERO
         ENDIF
      ENDDO

! Account for heat transfer between the discrete particles and the gas phase.
      IF(DES_CONTINUUM_COUPLED) CALL DES_Hgm(S_C, S_P)

! calculate the convection-diffusion terms
      CALL CONV_DIF_PHI (T_g, K_G, DISCRETIZE(6), U_G, V_G, W_G, &
         CpxFlux_E, CpxFlux_N, CpxFlux_T, 0, A_M, B_M, IER)

! calculate standard bc
      CALL BC_PHI (T_g, BC_T_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, 0, A_M, B_M, IER)

! set the source terms in a and b matrix equation form
      CALL SOURCE_PHI (S_P, S_C, EP_G, T_G, 0, A_M, B_M, IER)

! add point sources
      IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (T_g, PS_T_g, &
         PS_CpxMFLOW_g, 0, A_M, B_M, IER)

      DO M = 1, TMP_SMAX
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sE(IJK,M)
               ELSE   ! M=M_AM is the only phase for which virtual mass is added
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sSE(IJK)
               ENDIF
            ENDIF

            IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sN(IJK,M)
               ELSE
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sSN(IJK)
               ENDIF
            ENDIF

            IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sT(IJK,M)
               ELSE
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sST(IJK)
               ENDIF
            ENDIF

            IF (FLUID_AT(IJK)) THEN
               APO = ROP_SO(IJK,M)*C_PS(IJK,M)*VOL(IJK)*ODT
               S_P(IJK) = APO + S_RPS(IJK,M)*VOL(IJK)
               S_C(IJK) = APO*T_SO(IJK,M) - HOR_S(IJK,M)*VOL(IJK) + &
                  S_RCS(IJK,M)*VOL(IJK)
! Lane sub-grid models - 07/2015
               IF (SG_CYL_ENERGY) THEN
                  IF (EP_S(IJK,1) > DIL_EP_S .AND. EP_S(IJK,2) > ZERO)  THEN
                     ! Variables for Nusselt number calculation
                     EP_G_STAR = EP_G(IJK)/(1.0-EP_S(IJK,2))
                     EP_S_STAR = EP_S(IJK,1)/(1.0-EP_S(IJK,2))
                     V_T = D_P(IJK,1)**2*GRAVITY*(RO_S(IJK,1)-RO_G(IJK))/(18.0*MU_G(IJK))
                     L_STAR = V_T**2/GRAVITY
                     V_S_STAR = SQRT(U_S(IJK,1)**2 + V_S(IJK,1)**2 + W_S(IJK,1)**2)/V_T
                     PE_CS = RO_S(IJK,1)*C_PS(IJK,1)*V_T*L_STAR/K_S(IJK,1)
                     ! Nusselt number and heat transfer coefficient
                     NU_CS = EP_S_STAR**0.163 * V_S_STAR**0.0623 * (0.0777 + 0.448*(SG_CYL_D/SG_CYL_A)**3.55) * PE_CS**0.668
                     H_CS = NU_CS*K_S(IJK,1)/L_STAR
                     ! Area calculation and correction
                     L_R = (DX(I)*DY(J))/SG_CYL_A**2
                     A_C = 2.0*PI*(L_R*SG_CYL_D) ! 2D/3D Correction factor
                     ! Source term
                     T_SUSP = EP_G_STAR*T_G(IJK) + EP_S_STAR*T_S(IJK,1)
                     Q_CS = H_CS*A_C*(SG_CYL_T - T_SUSP)
                     S_C(IJK) = S_C(IJK) + Q_CS
                  ENDIF
               ENDIF
               VXGAMA(IJK,M) = GAMA_GS(IJK,M)*VOL(IJK)
               EPS(IJK) = EP_S(IJK,M)
               IF(USE_MMS) S_C(IJK) = S_C(IJK) + MMS_T_S_SRC(IJK)*VOL(IJK)
            ELSE
               S_P(IJK) = ZERO
               S_C(IJK) = ZERO
               VXGAMA(IJK,M) = ZERO
               EPS(IJK) = ZERO
               IF(USE_MMS) EPS(IJK) = EP_S(IJK,M)
            ENDIF
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
         CALL CONV_DIF_PHI (T_s(1,M), K_S(1,M), DISCRETIZE(6), &
            U_S(1,M), V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, &
            CpxFlux_T, M, A_M, B_M, IER)

! calculate standard bc
         CALL BC_PHI (T_s(1,M), BC_T_S(1,M), BC_TW_S(1,M), &
            BC_HW_T_S(1,M), BC_C_T_S(1,M), M, A_M, B_M, IER)

! set the source terms in a and b matrix equation form
         CALL SOURCE_PHI (S_P, S_C, EPS, T_S(1,M), M, A_M, B_M, IER)

! Add point sources.
         IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (T_s(:,M), PS_T_s(:,M),&
            PS_CpxMFLOW_s(:,M), M, A_M, B_M, IER)

      ENDDO   ! end do (m=1,tmp_smax)

! use partial elimination on interphase heat transfer term
      IF (TMP_SMAX > 0 .AND. .NOT.USE_MMS) &
        CALL PARTIAL_ELIM_S (T_G, T_S, VXGAMA, A_M, B_M, IER)

      CALL CALC_RESID_S (T_G, A_M, B_M, 0, NUM_RESID(RESID_T,0),&
         DEN_RESID(RESID_T,0), RESID(RESID_T,0), MAX_RESID(RESID_T,&
         0), IJK_RESID(RESID_T,0), ZERO, IER)

      CALL UNDER_RELAX_S (T_G, A_M, B_M, 0, UR_FAC(6), IER)

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      write(*,*) &
!         resid(resid_t, 0), max_resid(resid_t, 0), &
!         ijk_resid(resid_t, 0)


      DO M = 1, TMP_SMAX
         CALL CALC_RESID_S (T_S(1,M), A_M, B_M, M, NUM_RESID(RESID_T,M), &
            DEN_RESID(RESID_T,M), RESID(RESID_T,M), MAX_RESID(&
            RESID_T,M), IJK_RESID(RESID_T,M), ZERO, IER)

         CALL UNDER_RELAX_S (T_S(1,M), A_M, B_M, M, UR_FAC(6), IER)
      ENDDO

! set/adjust linear equation solver method and iterations
      CALL ADJUST_LEQ(RESID(RESID_T,0), LEQ_IT(6), LEQ_METHOD(6), &
         LEQI, LEQM, IER)
!      call test_lin_eq(a_m(1, -3, 0), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)

      CALL SOLVE_LIN_EQ ('T_g', 6, T_G, A_M, B_M, 0, LEQI, LEQM, &
         LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)
! Check for linear solver divergence.
      IF(ier == -2) Err_l(myPE) = 120

! bound temperature in any fluid or flow boundary cells
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.WALL_AT(IJK))&
            T_g(IJK) = MIN(TMAX, MAX(TMIN, T_g(IJK)))
      ENDDO

!      call out_array(T_g, 'T_g')

      DO M = 1, TMP_SMAX
         CALL ADJUST_LEQ (RESID(RESID_T,M), LEQ_IT(6), LEQ_METHOD(6), &
            LEQI, LEQM, IER)
!         call test_lin_eq(a_m(1, -3, M), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)
         CALL SOLVE_LIN_EQ ('T_s', 6, T_S(1,M), A_M, B_M, M, LEQI, &
            LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)

! Check for linear solver divergence.
         IF(ier == -2) Err_l(myPE) = 121

! bound temperature in any fluid or flow boundary cells
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.WALL_AT(IJK))&
               T_s(IJK, M) = MIN(TMAX, MAX(TMIN, T_s(IJK, M)))
         ENDDO
      ENDDO   ! end do (m=1, tmp_smax)

      call unlock_ambm
      call unlock_tmp_array
      call unlock_tmp_array1

! If the linear solver diverged, temperatures may take on unphysical
! values. To prevent them from propogating through the domain or
! causing failure in other routines, force an exit from iterate and
! reduce the time step.
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)


      RETURN
      END SUBROUTINE SOLVE_ENERGY_EQ
