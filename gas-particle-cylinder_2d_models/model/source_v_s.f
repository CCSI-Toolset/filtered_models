!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: SOURCE_V_s                                         C
!  Purpose: Determine source terms for V_s momentum eq. The terms      C
!     appear in the center coefficient and RHS vector.    The center   C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_V_S(A_M, B_M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE vshear
      USE compar
      USE sendrecv
      use kintheory
      USE ghdtheory
      USE drag
      USE cutcell
      USE quadric
      USE mms
      USE bodyforce
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IJKM, IJKN, &
                 IPJK, IJPK, IJKP, IMJPK, IJPKM
! Phase index
      INTEGER :: M, MM, L
! Internal surface
      INTEGER :: ISV
! Pressure at north cell
      DOUBLE PRECISION :: PgN
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp, epse, epsw, epsn, epss, &
                          epst, epsb, epsMix, epsMixN
      DOUBLE PRECISION :: SUM_EPS_CP
! Average density
      DOUBLE PRECISION :: ROPSA
! Average density difference
      DOUBLE PRECISION :: dro1, dro2, droa
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp, Sdps
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vmt, Vbf, Vmttmp
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, Vgn, Vgs, Uge, Ugw, Vge,&
                          Vgw, Wgt, Wgb, Vgt, Vgb
      DOUBLE PRECISION :: F_vir
! terms for shear (loezos)
      DOUBLE PRECISION :: VSH_n,VSH_s,VSH_e,VSH_w,VSH_p,Source_conv
      DOUBLE PRECISION :: SRT
! Sarkar sub-grid models - 07/2013
      DOUBLE PRECISION dragForceCylSuspY
      DOUBLE PRECISION solidsRatio
      DOUBLE PRECISION betaStarCylSuspY, gammaStarCylSuspY, vSolidStar
      DOUBLE PRECISION L_star, D_cyl_star, a_cyl_star
!-----------------------------------------------


      DO M = 1, MMAX
        IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
           (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

          IF (MOMENTUM_Y_EQ(M)) THEN


!$omp  parallel do default(shared)                                   &
!$omp  private( I, J, K, IJK, IJKN, IMJK, IPJK, IJMK, IJPK, IMJPK,   &
!$omp           IJKM, IJPKM, IJKP, ISV, epsMix, epsMixN, EPStmp,     &
!$omp           EPSA, EPSw, EPSe, EPSn, EPSs, EPSt, EPSb, SUM_EPS_CP,&
!$omp           PGN, Sdp, Sdps, MM, L, ROPSA, V0, ROP_MA,            &
!$omp           Vgn, Vgs, Uge, Ugw, Vge, Vgw, Wgt, Wgb, Vgt, Vgb,    &
!$omp           F_vir, VMT, VMTtmp, DRO1, DRO2, DROA, Vbf,           &
!$omp           avgRop, Ghd_drag, HYS_drag, avgDrag,                 &
!$omp           VSH_n, VSH_s, VSH_e, VSH_w, VSH_p, Source_conv, SRT)
            DO IJK = ijkstart3, ijkend3

! Skip walls where some values are undefined.
                IF(WALL_AT(IJK)) cycle

                I = I_OF(IJK)
                J = J_OF(IJK)
                K = K_OF(IJK)
                IMJK = IM_OF(IJK)
                IPJK = IP_OF(IJK)
                IJMK = JM_OF(IJK)
                IJPK = JP_OF(IJK)
                IMJPK = IM_OF(IJPK)
                IJKM = KM_OF(IJK)
                IJKN = NORTH_OF(IJK)
                IJPKM = KM_OF(IJPK)
                IJKP = KP_OF(IJK)

                IF (KT_TYPE_ENUM == GHD_2007) THEN
                  EPStmp = ZERO
                  epsMix = ZERO
                  epsMixN= ZERO
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J)
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixN to be used for model B
                    epsMixN = epsMixN + EP_S(IJKN,L)
                    IF(IP_AT_N(IJK)) THEN
                    V_S(IJK,L) = ZERO
                    ELSEIF(SIP_AT_N(IJK)) THEN
                       ISV = IS_ID_AT_N(IJK)
                       V_S(IJK,L) = IS_VEL_S(ISV,L)
                    ENDIF
                  ENDDO
                  EPSA = EPStmp
                ELSE
                  EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
                ENDIF

! Impermeable internal surface
                IF (IP_AT_N(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO

! Semi-permeable internal surface
                ELSEIF (SIP_AT_N(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  ISV = IS_ID_AT_N(IJK)
                  B_M(IJK,M) = -IS_VEL_S(ISV,M)

! Dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                      EPSw = ZERO
                      EPSe = ZERO
                      EPSn = ZERO
                      EPSs = ZERO
                      EPSt = ZERO
                      EPSb = ZERO
                      DO L = 1, SMAX
                        EPSw = EPSw + EP_S(WEST_OF(IJK),L)
                        EPSe = EPSe + EP_S(EAST_OF(IJK),L)
                        EPSn = EPSn + EP_S(NORTH_OF(IJK),L)
                        EPSs = EPSs + EP_S(SOUTH_OF(IJK),L)
                        IF(.NOT. NO_K) THEN
                          EPSt = EPSt + EP_S(TOP_OF(IJK),L)
                          EPSb = EPSb + EP_S(BOTTOM_OF(IJK),L)
                        ENDIF
                      ENDDO
                  ELSE
                      EPSw = EP_S(WEST_OF(IJK),M)
                      EPSe = EP_S(EAST_OF(IJK),M)
                      EPSn = EP_S(NORTH_OF(IJK),M)
                      EPSs = EP_S(SOUTH_OF(IJK),M)
                      IF(.NOT. NO_K) THEN
                        EPSt = EP_S(TOP_OF(IJK),M)
                        EPSb = EP_S(BOTTOM_OF(IJK),M)
                      ENDIF
                  ENDIF
! using the average boundary cell values to compute V_s (sof, Aug 23 2005)
                  IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE
                  IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE
                  IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE
                  IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF(.NOT. NO_K) THEN
                    IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE
                    IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE
                  ENDIF
                  IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
                    A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
                    B_M(IJK,M) = -V_S(IJK,M)
                  ELSE
                    A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
                  ENDIF

! Cartesian grid implementation
               ELSEIF (BLOCKED_V_CELL_AT(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO

! Normal case
               ELSE

! Surface forces:
! Pressure term
                  PGN = P_G(IJKN)
                  IF (CYCLIC_Y_PD) THEN
! CYCLIC_AT_N Flag is not set correctly in DMP and causes issues. This
! is avoided by using the DMP cyclic map. The flags need fixed.
!                     IF (CYCLIC_AT_N(IJK)) PGN = P_G(IJKN) - DELP_Y
                     IF (JMAP(J_OF(IJK)).EQ.JMAX1)PGN = P_G(IJKN) - DELP_Y
                  ENDIF
                  IF (MODEL_B) THEN
                     SDP = ZERO
                  ELSE
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        SDP = -P_SCALE*EPSA*(PGN - P_G(IJK))*AXZ(IJK)
                     ELSE
                        SDP = -P_SCALE*EPSA*(PGN * A_VPG_N(IJK)  - P_G(IJK) * A_VPG_S(IJK) )
                     ENDIF
                  ENDIF

                  IF (CLOSE_PACKED(M)) THEN
                     IF(SMAX > 1 .AND. KT_TYPE_ENUM /= GHD_2007) THEN
                        SUM_EPS_CP=0.0
                        DO MM=1,SMAX
                          IF (CLOSE_PACKED(MM))&
                            SUM_EPS_CP=SUM_EPS_CP+AVG_Y(EP_S(IJK,MM),EP_S(IJKN,MM),J)
                        ENDDO
                        SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                        SDPS = - ((P_S(IJKN,M)-P_S(IJK,M))+(EPSA/SUM_EPS_CP)* &
                           (P_STAR(IJKN)-P_STAR(IJK)))*AXZ(IJK)
                     ELSE
                        IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                          SDPS = - ((P_S(IJKN,M)-P_S(IJK,M))+(P_STAR(IJKN)-P_STAR(IJK)))*AXZ(IJK)
                        ELSE
                          SDPS = - ((P_S(IJKN,M)* A_VPG_N(IJK)-P_S(IJK,M)* A_VPG_S(IJK)) &
                               +(P_STAR(IJKN)* A_VPG_N(IJK)-P_STAR(IJK)* A_VPG_S(IJK)))
                        ENDIF
                     ENDIF
                  ELSE
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        SDPS = -(P_S(IJKN,M)-P_S(IJK,M))*AXZ(IJK)
                     ELSE
                        SDPS = -(P_S(IJKN,M) * A_VPG_N(IJK)-P_S(IJK,M) * A_VPG_S(IJK))
                     ENDIF
                  ENDIF


                  IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
! Volumetric forces
                     ROPSA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
! Previous time step
                     V0 = AVG_Y(ROP_SO(IJK,M),ROP_SO(IJKN,M),J)*ODT
! Added mass implicit transient term {Cv eps rop_g dV/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
                       ROP_MA = AVG_Y(ROP_g(IJK)*EP_s(IJK,M),ROP_g(IJKN)*EP_s(IJKN,M),J)
                       V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ELSE
! Volumetric forces
                     ROPSA =  (VOL(IJK)*ROP_S(IJK,M) +&
                        VOL(IJKN)*ROP_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
! Previous time step
                     V0 = (VOL(IJK)*ROP_SO(IJK,M) + &
                        VOL(IJKN)*ROP_SO(IJKN,M))*ODT/&
                        (VOL(IJK) + VOL(IJKN))
! Added mass implicit transient term {Cv eps rop_g dV/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
                        ROP_MA =  (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + &
                           VOL(IJKN)*ROP_g(IJKN)*EP_s(IJKN,M))/&
                           (VOL(IJK) + VOL(IJKN))
                        V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dvg/dt to virtual mass term
                  F_vir = ZERO
                  IF(Added_Mass .AND. M==M_AM .AND.&
                     (.NOT.CUT_V_TREATMENT_AT(IJK))) THEN
                     F_vir = ( (V_g(IJK) - V_gO(IJK)) )*ODT*VOL_V(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
                     Vgs = AVG_Y_N(V_G(IJMK),V_G(IJK))
                     Vgn = AVG_Y_N(V_G(IJK),V_G(IJPK))
                     Uge = AVG_Y(U_G(IJK),U_G(IJPK),J)
                     Ugw = AVG_Y(U_G(IMJK),U_G(IMJPK),J)
                     Vge = AVG_X(V_G(IJK),V_G(IPJK),IP1(I))
                     Vgw = AVG_X(V_G(IMJK),V_G(IJK),I)
                     IF(DO_K) THEN
                        Wgt = AVG_Y(W_g(IJK),W_g(IJPK),J)
                        Wgb = AVG_Y(W_g(IJKM),W_g(IJPKM),J)
                        Vgt = AVG_Z(V_g(IJK),V_g(IJKP),KP1(K))
                        Vgb = AVG_Z(V_g(IJKM),V_g(IJK),K)
                        F_vir = F_vir + AVG_Z_T(Wgb,Wgt)*&
                           OX(I) * (Vgt - Vgb)*AXY(IJK)
                     ENDIF
! adding convective terms (U dV/dx + V dV/dy) to virtual mass; W dV/dz added above.
                     F_vir = F_vir + V_g(IJK)*(Vgn - Vgs)*AXZ(IJK) + &
                        AVG_X_E(Ugw,Uge,IP1(I))*(Vge - Vgw)*AYZ(IJK)
                     F_vir = F_vir * Cv * ROP_MA
                  ENDIF

! Interphase mass transfer
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                     VMTtmp = ZERO
                     DO L = 1,SMAX
                        VMTtmp = VMTtmp + AVG_Y(SUM_R_S(IJK,L),SUM_R_S(IJKN,L),J)
                     ENDDO
                     VMT = VMTtmp
                  ELSE
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        VMT = AVG_Y(SUM_R_S(IJK,M),SUM_R_S(IJKN,M),J)
                     ELSE
                        VMT = (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IJKN)*&
                           SUM_R_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                     ENDIF
                  ENDIF

! Body force
                  IF (MODEL_B) THEN
                     IF (KT_TYPE_ENUM == GHD_2007) THEN
                       DRO1 = ROP_S(IJK,M)  - RO_G(IJK) *epsMix
                       DRO2 = ROP_S(IJKN,M) - RO_G(IJKN)*epsMixN
                       DROA = AVG_Y(DRO1,DRO2,J)
                       VBF = DROA*BFY_S(IJK,M)
                     ELSE
                       DRO1 = (RO_S(IJK,M)-RO_G(IJK))*EP_S(IJK,M)
                       DRO2 = (RO_S(IJK,M)-RO_G(IJKN))*EP_S(IJKN,M)
                       DROA = AVG_Y(DRO1,DRO2,J)
                       VBF = DROA*BFY_S(IJK,M)
                     ENDIF
                  ELSE
                     VBF = ROPSA*BFY_S(IJK,M)
                  ENDIF

! Additional force for GHD from drag force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                    DO L = 1,SMAX
                      avgRop = AVG_Y(ROP_S(IJK,L),ROP_S(IJKN,L),J)
                      if(avgRop > ZERO) Ghd_drag = Ghd_drag -&
                           AVG_Y(F_GS(IJK,L),F_GS(IJKN,L),J) * &
                           JoiY(IJK,L) / avgRop
                    ENDDO
                  ENDIF

! Additional force for HYS drag force, do not use with mixture GHD theory
                  HYS_drag = ZERO
                  avgDrag = ZERO
                  IF (DRAG_TYPE_ENUM .EQ. HYS .AND. KT_TYPE_ENUM /= GHD_2007) THEN
                     DO L = 1,MMAX
                        IF (L /= M) THEN
                           avgDrag = AVG_Y(beta_ij(IJK,M,L),beta_ij(IJKN,M,L),J)
                           HYS_drag = HYS_drag - avgDrag * (V_g(ijk) - V_s(IJK,L))
                        ENDIF
                     ENDDO
                  ENDIF

! Source terms from convective mom. flux (loezos)
                  IF (SHEAR) THEN
                    SRT=(2d0*V_sh/XLENGTH)
                    VSH_p=VSH(IJK)
                    VSH_n=VSH_p
                    VSH_s=VSH_p
                    VSH_e=VSH(IJK)+SRT*1d0/oDX_E(I)
                    VSH_w=VSH(IJK)-SRT*1d0/oDX_E(IM1(I))
                    Source_conv=A_M(IJK,N,m)*VSH_n+A_M(IJK,S,m)*VSH_s&
                      +A_M(IJK,W,m)*VSH_w+A_M(IJK,E,m)*VSH_e&
                      -(A_M(IJK,N,m)+A_M(IJK,S,m)+A_M(IJK,W,m)+A_M(IJK,E,m))&
                      *VSH_p
                  ELSE
                    Source_conv=0d0
                  ENDIF

! Sarkar sub-grid model - 07/2015
                  IF (SG_CYL_HYDRO) THEN
                     dragForceCylSuspY = ZERO
                     IF (M == 1 .AND. EP_S(IJK,2) > ZERO .AND. FLUID_AT(IJK)) THEN !drag exerted only on the flowing solids phase
                        V_t = GRAVITY*D_p0(1)*D_p0(1)*(RO_S(IJK,1)-RO_g(IJK))/(18.0d0*MU_G(IJK))
                        L_star = V_t*V_t/GRAVITY
                        D_cyl_star = SG_CYL_D/L_star
                        a_cyl_star = SG_CYL_a/L_star
                        
                        B1 = ( - 0.1106*(D_cyl_star**4) + 1.047*(D_cyl_star**3) &
                             - 2.354*(D_cyl_star**2) + 1.957*D_cyl_star ) & 
                             / ( a_cyl_star**2 - 22.74*a_cyl_star + 134.0 )
                        B2 = ( - 6.273*(D_cyl_star**3) + 40.86*(D_cyl_star**2) ) &
                             / ( a_cyl_star**2 - 26.86*a_cyl_star + 196.3 )
                        
                        solidsRatio = EP_S(IJK,1)/(ONE - EP_S(IJK,2)) !used repeatedly

                        betaStarCylSuspY = B1*(solidsRatio**2) &
                             / (ONE + B2*(solidsRatio**2))
                        gammaStarCylSuspY = (EP_S(IJK,2)/(ONE - EP_S(IJK,2))) &
                             * G1*solidsRatio/(ONE + G2*solidsRatio)

                        ! ijmk = JM_OF(ijk) !already calculated
                        vSolidStar = HALF*( V_s(ijk,M) + V_s(IJMK,M))/V_t
                        ! dimensionless (avge) cell center solids Y-velocity 
                        
                        ! filtered drag force due to cylinders on solids, per unit volume
                        dragForceCylSuspY = RO_S(IJK,1)*GRAVITY*(ONE - EP_S(IJK,2)) & 
                             *(-betaStarCylSuspY*vSolidStar*abs(vSolidStar) - gammaStarCylSuspY )

                     ENDIF
                  ENDIF
! /Sarkar sub-grid model


! Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
                     A_M(IJK,N,M)+A_M(IJK,S,M)+A_M(IJK,T,M)+&
                     A_M(IJK,B,M)+(V0+ZMAX(VMT))*VOL_V(IJK))

                  B_M(IJK,M) = B_M(IJK,M) - (SDP + SDPS + &
                        TAU_V_S(IJK,M) + Source_conv + F_vir + &
                        ( (V0+ZMAX((-VMT)))*V_SO(IJK,M) + &
                        VBF + HYS_drag)*VOL_V(IJK) )
                        
! Sarkar sub-grid models - 07/2015
! Add the contridution due to the subgrid cylinder-solids drag
                  IF (SG_CYL_HYDRO) THEN
                     B_M(IJK,M) = B_M(IJK,M) - dragForceCylSuspY*VOL_V(IJK)
                  ENDIF
                  
! MMS Source term
                  IF(USE_MMS) B_M(IJK,M) = &
                     B_M(IJK,M) - MMS_V_S_SRC(IJK)*VOL_V(IJK)

                  IF (KT_TYPE_ENUM == IA_2005) THEN
                     B_M(IJK,M) = B_M(IJK,M) - KTMOM_V_S(IJK,M)
                  ELSEIF (KT_TYPE_ENUM == GHD_2007) THEN
                     B_M(IJK,M) = B_M(IJK,M) - Ghd_drag*VOL_V(IJK)
                  ENDIF

                ENDIF   ! end branching on cell type (sip/ip/dilute/block/else branches)
            ENDDO   ! end do loop over ijk
!$omp end parallel do

! modifications for cartesian grid implementation
            IF(CARTESIAN_GRID) CALL CG_SOURCE_V_S(A_M, B_M, M, IER)
! modifications for bc
            CALL SOURCE_V_S_BC (A_M, B_M, M, IER)
            IF(CARTESIAN_GRID) CALL CG_SOURCE_V_S_BC(A_M, B_M, M, IER)

          ENDIF   ! end if (momentum_y_eq)
        ENDIF   ! end if for GHD Theory
      ENDDO   ! end do loop over mmax

      RETURN
      END SUBROUTINE SOURCE_V_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_V_s_BC                                           C
!  Purpose: Determine source terms for V_s momentum eq. The terms      C
!     appear in the center coefficient and RHS vector.    The center   C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Comments: see source_v_g_bc for more detailed in-code comments      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_V_S_BC(A_M, B_M, M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER ::  L
! Indices
      INTEGER ::  I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                  IM, KM, IJKS, IJMK, IJPK
!-----------------------------------------------

! Setting the default boundary conditions
      IF (DO_K) THEN
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = -ONE
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ONE
                  A_M(IJK,B,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO

         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = -ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSEIF (FS_WALL_AT(IJK)) THEN
                  A_M(IJK,E,M) = ZERO
                  A_M(IJK,W,M) = ZERO
                  A_M(IJK,N,M) = ZERO
                  A_M(IJK,S,M) = ZERO
                  A_M(IJK,T,M) = ZERO
                  A_M(IJK,B,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = -ONE
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ONE
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      I1 = IMAX2
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = -ONE
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ONE
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO
! End setting the default boundary conditions

! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ZERO
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              A_M(IJK,E,M) = -ONE
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              A_M(IJK,W,M) = -ONE
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                              A_M(IJK,T,M) = -ONE
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                              A_M(IJK,B,M) = -ONE
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

            ELSEIF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ZERO
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              A_M(IJK,E,M) = ONE
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              A_M(IJK,W,M) = ONE
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                              A_M(IJK,T,M) = ONE
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                              A_M(IJK,B,M) = ONE
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           IM = IM1(I)
                           KM = KM1(K)
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ZERO
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,E,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_VW_S(L,M)
                              ELSE
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(I))
                                 A_M(IJK,E,M) = -(HALF*BC_HW_S(L,M)-ODX_E(I))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,W,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_VW_S(L,M)
                              ELSE
                                 A_M(IJK,W,M) = -(HALF*BC_HW_S(L,M)-ODX_E(IM))
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(IM))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,T,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_VW_S(L,M)
                              ELSE
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*OX&
                                    (I))
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*OX&
                                    (I))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,B,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_VW_S(L,M)
                              ELSE
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX(I))
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX(I))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M)
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

! Setting flow boundary conditions
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'S') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ONE
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF

            ELSEIF (BC_TYPE(L) == 'OUTFLOW') THEN
               IF (BC_PLANE(L) == 'S') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,E,M) = ZERO
                           A_M(IJK,W,M) = ZERO
                           A_M(IJK,N,M) = ZERO
                           A_M(IJK,S,M) = ONE
                           A_M(IJK,T,M) = ZERO
                           A_M(IJK,B,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IJMK = JM_OF(IJK)
                           A_M(IJMK,E,M) = ZERO
                           A_M(IJMK,W,M) = ZERO
                           A_M(IJMK,N,M) = ZERO
                           A_M(IJMK,S,M) = ONE
                           A_M(IJMK,T,M) = ZERO
                           A_M(IJMK,B,M) = ZERO
                           A_M(IJMK,0,M) = -ONE
                           B_M(IJMK,M) = ZERO
                        END DO
                     END DO
                  END DO
               ELSEIF (BC_PLANE(L) == 'N') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IJPK = JP_OF(IJK)
                           A_M(IJPK,E,M) = ZERO
                           A_M(IJPK,W,M) = ZERO
                           A_M(IJPK,N,M) = ONE
                           A_M(IJPK,S,M) = ZERO
                           A_M(IJPK,T,M) = ZERO
                           A_M(IJPK,B,M) = ZERO
                           A_M(IJPK,0,M) = -ONE
                           B_M(IJPK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF

! Setting bc that are not ns, fs, psw, p_inflow, p_outflow, or outflow
            ELSE
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        A_M(IJK,E,M) = ZERO
                        A_M(IJK,W,M) = ZERO
                        A_M(IJK,N,M) = ZERO
                        A_M(IJK,S,M) = ZERO
                        A_M(IJK,T,M) = ZERO
                        A_M(IJK,B,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = -V_S(IJK,M)
                        IF (BC_PLANE(L) == 'S') THEN
                           IJKS = SOUTH_OF(IJK)
                           A_M(IJKS,E,M) = ZERO
                           A_M(IJKS,W,M) = ZERO
                           A_M(IJKS,N,M) = ZERO
                           A_M(IJKS,S,M) = ZERO
                           A_M(IJKS,T,M) = ZERO
                           A_M(IJKS,B,M) = ZERO
                           A_M(IJKS,0,M) = -ONE
                           B_M(IJKS,M) = -V_S(IJKS,M)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF   ! end if (bc_type)
         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN
      END SUBROUTINE SOURCE_V_S_BC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_V_s                                              C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal, A. Srivastava,                 Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE JJ_BC_V_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE output
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Boundary condition
      INTEGER, INTENT(IN) :: L
! Indices
      INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2
! Solids phase index
      INTEGER, INTENT(IN) ::  M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IM, KM, IJPK
! coefficients for granular bc
      DOUBLE PRECISION :: hw, gw, cw
!-----------------------------------------------

      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

               IJK = FUNIJK(I,J,K)
               IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
               IM = IM1(I)
               KM = KM1(K)
               A_M(IJK,E,M) = ZERO
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO
               A_M(IJK,S,M) = ZERO
               A_M(IJK,T,M) = ZERO
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IF (FLUID_AT(EAST_OF(IJK))) THEN
                  IJPK = JP_OF(EAST_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(EAST_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,E,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
! Setting the wall velocity based on Johnson and Jackson B.C which are
! also modified to include frictional effects if requested
                        CALL CALC_GRBDRY (IJK, EAST_OF(IJK), 'E', 'V', &
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
! Setting the wall velocity equal to the user specified value
! (if wall velocity = 0 then it is equivalent to no slip wall)
                        GW = 0D0
                        HW = 1D0
                        CW = BC_VW_S(L,M)
                     ELSE
! Setting the wall velocity equal to the adjacent fluid cell
! velocity (i.e. zero flux/free slip)
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,E,M) = -(HALF*HW - ODX_E(I)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODX_E(I)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                  IJPK = JP_OF(WEST_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(WEST_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,W,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, WEST_OF(IJK), 'W', 'V',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_VW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,W,M) = -(HALF*HW - ODX_E(IM)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODX_E(IM)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                  IJPK = JP_OF(TOP_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,T,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'V',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_VW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX(I)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX(I)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                  IJPK = JP_OF(BOTTOM_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,B,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'V',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_VW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX(I)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX(I)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE JJ_BC_V_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_V_S                                        C
!  Purpose: Adds point sources to the solids V-Momentum equations.     C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_V_S(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use constant
      use fldvar
      use geometry
      use indices
      use physprop
      use ps
      use run
      use functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lJN, lJS
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

      do M=1, MMAX

      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_V_s(PSV,M)) < small_number) cycle PS_LP

         if(PS_V_s(PSV,M) < 0.0d0) then
            lJS = PS_J_S(PSV) - 1
            lJN = PS_J_N(PSV) - 1
         else
            lJS = PS_J_S(PSV)
            lJN = PS_J_N(PSV)
         endif

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = lJS, lJN
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle

            if(A_M(IJK,0,M) == -ONE .AND.                              &
               B_M(IJK,M) == -V_s(IJK,M)) then
               B_M(IJK,M) = -PS_V_s(PSV,M) * PS_VEL_MAG_S(PSV,M)
            else
               pSource = PS_MASSFLOW_S(PSV,M) *                        &
                  (VOL(IJK)/PS_VOLUME(PSV))

               B_M(IJK,M) = B_M(IJK,M) - pSource *                     &
                  PS_V_s(PSV,M) * PS_VEL_MAG_S(PSV,M)
            endif

         enddo
         enddo
         enddo


         enddo PS_LP

      enddo ! do M=1, MMAX

      RETURN
      END SUBROUTINE POINT_SOURCE_V_S
