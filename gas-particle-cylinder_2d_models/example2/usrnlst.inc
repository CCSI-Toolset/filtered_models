!
!  First add the variable name in the NAMELIST section.
!  Second define the user-defined namelist variables (e.g., DOUBLE PRECISION, 
!  INTEGER, LOGICAL, or CHARACTER.) in usr_mod.f. The user need not include this
!  file anywhere; it is already included in read_namelist.

!  A dummy variable is defined, so that the code can be compiled with the
!  template include file.  
      NAMELIST / USR_INPUT_DATA / DUMMY_DP, USE_CONC, REAC_E, REAC_DH, REAC_DS, REAC_LOGZETA
