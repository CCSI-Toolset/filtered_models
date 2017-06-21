!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_ss                                                 C
!  Purpose: This module computes the coefficient of drag between       C
!     two solids phases (M and L).                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     Gera, D., Syamlal, M., O'Brien T.J. 2004. International Journal  C
!        of Multiphase Flow, 30, p419-428.                             C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DRAG_SS(L, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE discretelement
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE rdf
      USE sendrecv
      USE run  ! Lane/Sarkar sub-grid models - 07/2015
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Index of solids phases
      INTEGER, INTENT(IN) :: L, M
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! indices
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
      INTEGER :: CM, DM
! index for storing solids-solids drag coefficients in the upper
! triangle of the matrix
      INTEGER :: LM
! cell center value of U_sm, U_sl, V_sm, V_sl, W_sm, W_sl
      DOUBLE PRECISION :: USCM, USCL, VSCM, VSCL, WSCM, WSCL
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VREL
! particle diameters of phase M and phase L
      DOUBLE PRECISION :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION :: RO_M, RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION :: G0_ML
! void fraction and solids volume fraction
      DOUBLE PRECISION :: EPg, EPs
! sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
!-----------------------------------------------

      LM = FUNLM(L,M)

      DO IJK = ijkstart3, ijkend3

         IF (.NOT.WALL_AT(IJK)) THEN
! Evaluate at all flow boundaries and fluid cells
! This is unlike the fluid-solid drag coefficient, which is only
! evluated in fluid cells and pressure inflow cells

            I = I_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! calculating velocity components at i, j, k (cell center)
            USCL = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I)
            VSCL = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L))
            WSCL = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L))

            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! magnitude of solids-solids relative velocity
            VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + &
                        (WSCL - WSCM)**2)

! setting aliases for easy reference
            D_PM = D_P(IJK,M)
            D_PL = D_P(IJK,L)
            RO_M = RO_S(IJK,M)
            RO_L = RO_S(IJK,L)

            IF (DES_CONTINUUM_HYBRID) THEN
! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
               EPSoDP = ZERO
               DO CM = 1, MMAX
                  EPS = EP_s(IJK, CM)
                  EPSoDP = EPSoDP + EPS / D_p(IJK,CM)
               ENDDO
               DO DM = 1, DES_MMAX
                  EPS = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
                  EPSoDP = EPSoDP + EPS / DES_D_p0(DM)
               ENDDO
               EPg = EP_g(IJK)
               G0_ML = ONE/EPg + 3.0d0*EPSoDP*D_pM*D_PL / &
                  (EPg*EPg *(D_pM + D_pL))
            ELSE
               G0_ML = G_0(IJK,L,M)
            ENDIF

! determining the solids-solids 'drag coefficient'
            CALL DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L,G0_ML,VREL)

            F_SS(IJK,LM) = lDss*ROP_S(IJK,M)*ROP_S(IJK,L)

! Gera: accounting for particle-particle drag due to enduring contact in a
! close-packed system. Note the check for mmax >= 2 below is unnecessary
! since this routine will not be entered unless mmax >=2
            IF(CLOSE_PACKED(M) .AND. CLOSE_PACKED(L) .AND. &
              (MMAX >= 2))  F_SS(IJK,LM) = F_SS(IJK,LM) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

          ELSE   ! elseif (.not.wall_at(ijk))

            F_SS(IJK,LM) = ZERO

         ENDIF   ! end if (.not.wall_at(ijk))

! Sarkar subgrid models - 07/2015
         IF (SG_CYL_HYDRO) THEN
            IF (L == 2 .OR. M == 2) THEN
               ! If either solid phases represent the cylinders, then set the drag
               ! coeff to zero.  The cyl-soilds (aka cyl-suspension) drag is modeled
               ! as source terms in source_v_s.f & source_u_s.f
               F_SS(IJK,LM) = ZERO
            ENDIF
         ENDIF
! /Sarkar subgrid models

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE DRAG_SS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SS_SYAM                                            C
!  Purpose: Calculate the solids-solids drag coefficient between a     C
!           continuous solids phase and discrete solids                C
!                                                                      C
!  Literature/Document References:                                     C
!     M. Syamlal. 1987. The particle-particle drag term in a           C
!        multiparticle model of fluidization. Technical Report.        C
!        DOE/MC/21353-2373. Office of Fossil Energy, Morgantown        C
!        Energy Technology Center, Morgantown, West Virginia.          C
!                                                                      C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SS_SYAM(lDss,D_PM,D_PL,RO_M,RO_L, G0_ML, VREL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      use run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: ldss
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: D_PM
! particle diameter of solids phase L
      DOUBLE PRECISION, INTENT(IN) :: D_PL
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: RO_M
! particle density of solids phase L
      DOUBLE PRECISION, INTENT(IN) :: RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION, INTENT(IN) :: G0_ML
! relative velocity between solids phase m and l
      DOUBLE PRECISION, INTENT(IN) :: VREL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Sum of particle diameters
      DOUBLE PRECISION :: DPSUM
! Intermediate calculation
      DOUBLE PRECISION :: const
!-----------------------------------------------
! External functions
!-----------------------------------------------
!-----------------------------------------------

      DPSUM = D_PL + D_PM

      const = 3.d0*(ONE + C_E)*(PI/2.d0 + C_F*PI*PI/8.d0)*&
         DPSUM**2/(2.d0*PI*(RO_L*D_PL**3+RO_M*D_PM**3))

      ldss = const * G0_ML * VREL

      RETURN
      END SUBROUTINE DRAG_SS_SYAM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_ss_IA                                              !
!  Purpose: Compute the coefficient of drag between solids phase m     !
!           and solids phase l using Iddir Arastoopour (2005) kinetic  !
!           theory model                                               !
!                                                                      !
!  Literature/Document References:                                     !
!    Iddir, Y.H., "Modeling of the multiphase mixture of particles     !
!       using the kinetic theory approach," PhD Thesis, Illinois       !
!       Institute of Technology, Chicago, Illinois, 2004               !
!    Iddir, Y.H., & H. Arastoopour, "Modeling of Multitype particle    !
!      flow using the kinetic theory approach," AIChE J., Vol 51,      !
!      no. 6, June 2005                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DRAG_SS_IA(L, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE drag
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE kintheory
      USE param1
      USE physprop
      USE rdf
      USE sendrecv

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M, L
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
! Index for storing solids-solids drag coefficients
! in the upper triangle of the matrix
      INTEGER :: LM
! Particle diameters
      DOUBLE PRECISION :: D_PM, D_PL
! Sum of particle diameters
      DOUBLE PRECISION :: DPSUM
!
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, DPSUMo2, NU_PL, NU_PM
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R0p_lm, R2p_lm, R3p_lm, &
                          R4p_lm, R10p_lm, Bp_lm
      DOUBLE PRECISION :: Fss_ip, Fnus_ip, FTsM_ip, FTsL_ip, &
                          F_common_term
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3

          IF (.NOT.WALL_AT(IJK)) THEN

               LM = FUNLM(L,M)

               IF (M == L) THEN
                    F_SS(IJK,LM) = ZERO
                    Fnu_s_ip(IJK,M,L) = ZERO
                    FT_sM_ip(IJK,M,L) = ZERO
                    FT_sL_ip(IJK,M,L) = ZERO

               ELSE
                    D_PM = D_P(IJK,M)
                    D_PL = D_P(IJK,L)
                    DPSUM = D_PL + D_PM
                    M_PM = (Pi/6.d0) * D_PM**3 *RO_S(IJK,M)
                    M_PL = (Pi/6.d0) * D_PL**3 *RO_S(IJK,L)
                    MPSUM = M_PM + M_PL
                    DPSUMo2 = DPSUM/2.d0
                    NU_PM = ROP_S(IJK,M)/M_PM
                    NU_PL = ROP_S(IJK,L)/M_PL

                 IF(Theta_m(IJK,M) > ZERO .AND. Theta_m(IJK,L) > ZERO) THEN

                    Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/&
                          2.d0
                    Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-Theta_m(IJK,M) ))/&
                         (2.d0*MPSUM)
                    Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+M_PL*Theta_m(IJK,L) ))/&
                         (2.d0*MPSUM*MPSUM)

                    R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+ &
                              ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                              ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 * Dp_lm**4.5 ) )

                    R2p_lm = ( 1.d0/( 2.d0*Ap_lm**1.5 * Dp_lm*Dp_lm ) )+&
                              ( (3.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**3 ) )+&
                              ( (15.d0*Bp_lm**4)/( 2.d0*Ap_lm**3.5 * Dp_lm**4 ) )

                    R3p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) ) )+&
                              ( (21.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**2.5 * Dp_lm**4.5 ) )+&
                              ( (315.d0*Bp_lm**4)/( 8.d0 * Ap_lm**3.5 *Dp_lm**5.5 ) )

                    R4p_lm = ( 3.d0/( Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                              ( (35.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**3.5 * Dp_lm**4.5 ) )+&
                              ( (441.d0*Bp_lm**4)/( 8.d0 * Ap_lm**4.5 * Dp_lm**5.5 ) )

                    R10p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**2.5 ) )+&
                              ( (25.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**3.5 * Dp_lm**3.5 ) )+&
                              ( (1225.d0*Bp_lm**4)/( 24.d0* Ap_lm**4.5 * Dp_lm**4.5 ) )

                    F_common_term = (DPSUMo2*DPSUMo2/4.d0)*(M_PM*M_PL/MPSUM)*&
                         G_0(IJK,M,L)*(1.d0+C_E)*(M_PM*M_PL)**1.5

! Momentum source associated with relative velocity between solids
! phase m and solid phase l
                    Fss_ip = F_common_term*NU_PM*NU_PL*DSQRT(PI)*R2p_lm*&
                           (Theta_m(IJK,M)*Theta_m(IJK,L))**2

! Momentum source associated with the difference in the gradients in
! number density of solids phase m and all other solids phases
                    Fnus_ip = F_common_term*(PI*DPSUMo2/12.d0)*R0p_lm*&
                           (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

! Momentum source associated with the gradient in granular temperature
! of solid phase M
                    FTsM_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                              (Theta_m(IJK,M)**1.5 * Theta_m(IJK,L)**2.5) *&
                              (  (-1.5d0/12.d0*R0p_lm)+&
                              Theta_m(IJK,L)/16.d0*(  (-M_PM*R10p_lm) - &
                              ((5.d0*M_PL*M_PL*M_PM/(192.d0*MPSUM*MPSUM))*R3p_lm)+&
                              ((5.d0*M_PM*M_PL)/(96.d0*MPSUM)*R4p_lm*Bp_lm)  )  )

! Momentum source associated with the gradient in granular temperature
! of solid phase L ! no need to recompute (sof Aug 30 2006)
                    FTsL_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                             (Theta_m(IJK,L)**1.5 * Theta_m(IJK,M)**2.5) *&
                              (  (1.5d0/12.d0*R0p_lm)+&
                              Theta_m(IJK,M)/16.d0*(  (M_PL*R10p_lm)+&
                              (5.d0*M_PM*M_PM*M_PL/(192.d0*MPSUM*MPSUM)*R3p_lm)+&
                              (5.d0*M_PM*M_PL/(96.d0*MPSUM) *R4p_lm*Bp_lm)  )  )

                    F_SS(IJK,LM) = Fss_ip

                    Fnu_s_ip(IJK,M,L) = Fnus_ip

! WARNING: the following two terms have caused some convergence problems
! earlier. Set them to ZERO for debugging in case of converegence
! issues. (sof)
                    FT_sM_ip(IJK,M,L) = FTsM_ip  ! ZERO

                    FT_sL_ip(IJK,M,L) = FTsL_ip  ! ZERO
                 ELSE
                   F_SS(IJK,LM) = ZERO
                    Fnu_s_ip(IJK,M,L) = ZERO
                    FT_sM_ip(IJK,M,L) = ZERO
                    FT_sL_ip(IJK,M,L) = ZERO
                 ENDIF
               ENDIF
          ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DRAG_SS_IA

