      NAMELIST / INPUT_DATA /&

! 1.1 Run control section
! ---------------------------------------------------------------->>>
             RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME, TSTOP,&
	     Dt, DT_MAX, DT_MIN, DT_FAC, DETECT_STALL,&
             AUTO_RESTART,&	     
             ENERGY_EQ, SPECIES_EQ, CALL_USR, &
             MOMENTUM_X_EQ, MOMENTUM_Y_EQ, MOMENTUM_Z_EQ, &
	     Call_DQMOM, AGGREGATION_EFF,BREAKAGE_EFF, NScalar, &
	     Phase4Scalar, &
	     K_Epsilon, SIMONIN, AHMADI, JENKINS, & 
             FRICTION, SAVAGE, SCHAEFFER,&
             BLENDING_STRESS, TANH_BLEND, SIGM_BLEND, &
             SHEAR, V_sh,&
             DRAG_TYPE, MODEL_B, lam_HYS, drag_c1, drag_d1, &
             SUBGRID_TYPE, filter_size_ratio, SUBGRID_Wall, & 
             GRANULAR_ENERGY, KT_TYPE, RDF_TYPE, &
             YU_STANDISH, FEDORS_LANDEL, Added_Mass, M_AM, &

! see run_mod for details
             CHK_BATCHQ_END, BATCH_WALLCLOCK, TERM_BUFFER,&

! Lane/Sarkar subgrid models 2015
  	      SG_3D_HYDRO, SG_CYL_HYDRO, SG_CYL_ENERGY, &
	      SG_CYL_D, SG_CYL_A, SG_CYL_T, &

! ----------------------------------------------------------------<<<

! 1.2 Output control section
! ---------------------------------------------------------------->>>
             RES_DT, SPX_DT , OUT_DT, &
             USR_DT, USR_TYPE, USR_VAR, USR_FORMAT, USR_EXT,&
             USR_X_w , USR_X_e , USR_Y_s , USR_Y_n ,&
             USR_Z_b , USR_Z_t, NLOG, FULL_LOG, &
             RESID_STRING,GROUP_RESID, &
	     report_mass_balance_dt, nRR, &

! see output_mod for details
             DBGPRN_LAYOUT, ENABLE_DMP_LOG, &
! ----------------------------------------------------------------<<<


! 1.3 Physical and numerical parameters section
! ---------------------------------------------------------------->>>
!            (e, r_p, and C_e are same)
             GRAVITY, GRAVITY_X,GRAVITY_Y,GRAVITY_Z,P_ref, P_scale, &
             e, C_e, r_p, C_f, C, C_NAME, Phi, Phi_w, L_scale0, MU_gmax,&
             PHIP, e_w, EPS_f_min, V_ex, SEGREGATION_SLOPE_COEFFICIENT,&
	     ep_s_max, M_MAX,&
             MAX_NIT, NORM_g, NORM_s, &
	     TOL_RESID, TOL_RESID_T, TOL_RESID_X,&
	     TOL_RESID_Th, TOL_RESID_Scalar, TOL_DIVERGE, &
	     MAX_INLET_VEL_FAC, &
             DISCRETIZE, LEQ_IT, LEQ_METHOD, &
	     LEQ_SWEEP, LEQ_TOL, LEQ_PC, DO_TRANSPOSE, &
	     UR_FAC, UR_F_gs, UR_Kth_sml, &	     
             DEF_COR, FPFOI, C_FAC, &
             chi_scheme, CN_ON, &
             BC_JJ_M, PHIP_OUT_JJ, PHIP0, &
! ----------------------------------------------------------------<<<


! 1.4 Geometry and discretization section
! ---------------------------------------------------------------->>>
             COORDINATES, IMAX, DX, XMIN, XLENGTH, JMAX, DY, YLENGTH,&
             KMAX, DZ, ZLENGTH, NO_I, NO_J, NO_K,&
             CYCLIC_X, CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD, CYCLIC_Y_PD,&
             CYCLIC_Z_PD, &
! ----------------------------------------------------------------<<<

! 1.5 Gas phase section
! ---------------------------------------------------------------->>>
             RO_g0, MU_g0, MW_AVG, MW_g, K_g0, C_pg0, DIF_g0,&
! ----------------------------------------------------------------<<<

! 1.6 Solids phase section
! ---------------------------------------------------------------->>>
             MU_s0, D_p0, EP_star, MMAX, MW_s, CLOSE_PACKED, &
             K_s0, C_ps0, DIF_s0, NMAX_g, NMAX_s, &
             RO_s0, X_s0, RO_Xs0, INERT_SPECIES, &
             SOLIDS_MODEL, &
! ----------------------------------------------------------------<<<


! 1.7 Initial conditions section
! ---------------------------------------------------------------->>>
             IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t,&
             IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t,&
             IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s,&
             IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s,&
             IC_X_g, IC_X_s, IC_TYPE, IC_P_star, IC_L_scale,&
             IC_GAMA_Rg, IC_T_Rg, IC_GAMA_Rs, IC_T_Rs,&
             IC_Theta_m, IC_Scalar, IC_K_Turb_G, IC_E_Turb_G, &
	     IC_EP_s, IC_DES_FIT_TO_REGION, IC_PIC_CONST_NPC, &
	     IC_PIC_CONST_STATWT, &
! ----------------------------------------------------------------<<<

! 1.8 Boundary conditions section
! ---------------------------------------------------------------->>>
             BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t,&
             BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t,&
             BC_EP_g, BC_P_g, BC_RO_g, BC_ROP_g, BC_ROP_s,&
             BC_EP_s, BC_T_g, BC_T_s,&
             BC_U_g, BC_U_s, BC_V_g, BC_V_s, BC_W_g, BC_W_s, BC_VELMAG_G,BC_VELMAG_S,&
             BC_TYPE,BC_PO_APPLY_TO_DES, BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g, &
             BC_MI_AS_WALL_FOR_DES, &
             BC_MASSFLOW_s, BC_X_g, BC_X_s,&
             BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_l, &
             BC_Jet_gl,&
             DELP_X, DELP_Y, DELP_Z, U_g0, U_s0, V_g0, V_s0, W_g0, W_s0,&
             BC_hw_g, BC_hw_s,&
             BC_Uw_g, BC_Uw_s, BC_Vw_g, BC_Vw_s, BC_Ww_g, BC_Ww_s,&
             BC_hw_T_g, BC_hw_T_s, BC_Tw_g, BC_Tw_s, BC_C_T_g,&
             BC_C_T_s,&
             BC_hw_X_g, BC_hw_X_s, BC_Xw_g, BC_Xw_s, BC_C_X_g,&
             BC_C_X_s,&
             BC_JJ_PS,&
             BC_hw_Theta_m, BC_Thetaw_m, BC_C_Theta_m,&
             BC_Theta_m,&
             BC_hw_Scalar, BC_Scalarw, BC_C_Scalar,&
             BC_Scalar,&
	     BC_K_Turb_G, BC_E_Turb_G, &
	     Flux_g, &
             BC_PIC_MI_CONST_NPC,BC_PIC_MI_CONST_STATWT, &
! ----------------------------------------------------------------<<<

! 1.X Point Sources
! ---------------------------------------------------------------->>>
             PS_X_w, PS_X_e, PS_Y_s, PS_Y_n, PS_Z_b, PS_Z_t,&
             PS_I_w, PS_I_e, PS_J_s, PS_J_n, PS_K_b, PS_K_t,&
             PS_U_g, PS_V_g, PS_W_g, PS_MASSFLOW_g, &
             PS_U_s, PS_V_s, PS_W_s, PS_MASSFLOW_s, &
             PS_T_g, PS_T_s, PS_X_g, PS_X_s,&
! ----------------------------------------------------------------<<<

! 1.9 Internal surfaces section
! ---------------------------------------------------------------->>>
             IS_X_w, IS_X_e, IS_Y_s, IS_Y_n, IS_Z_b, IS_Z_t,&
             IS_I_w, IS_I_e, IS_J_s, IS_J_n, IS_K_b, IS_K_t,&
             IS_TYPE, IS_PC, IS_VEL_s,&
! ----------------------------------------------------------------<<<

! 1.10 Chemical reactions
! ---------------------------------------------------------------->>>
            USE_RRATES, &
            SPECIES_g, SPECIES_s, &
            SPECIES_ALIAS_g, SPECIES_ALIAS_s, &
            STIFF_CHEMISTRY, STIFF_CHEM_MAX_STEPS, &
! ----------------------------------------------------------------<<<

! 1.11 parallel related
! ---------------------------------------------------------------->>>
             NODESI , NODESJ , NODESK, &
             DEBUG_RESID, &
	     IS_SERIAL, USE_DOLOOP, minimize_dotproducts, &
	     icheck_bicgs, &
             solver_statistics, opt_parallel , &

! Additional Debug Flags
! ---------------------------------------------------------------->>>
             REPORT_NEG_DENSITY, &
! Distributed IO (version 1) related variables
             bDist_IO , bStart_with_one_RES , &
! ----------------------------------------------------------------<<<

! Legacy variables
! ---------------------------------------------------------------->>>
            NMAX, & ! Legacy access
            SPECIES_NAME, &  ! Legacy variable
! ----------------------------------------------------------------<<<


! Misc
! ---------------------------------------------------------------->>>
             bWrite_netCDF
! ----------------------------------------------------------------<<<

