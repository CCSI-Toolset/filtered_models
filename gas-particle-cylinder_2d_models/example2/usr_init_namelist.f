!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_INIT_NAMELIST                                      C
!  Purpose: initialize user_defined NAMELIST variables                 C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR_INIT_NAMELIST 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      Use usr
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      USE_CONC = .FALSE.
!   1. Dry CO2 Adsorption
!      2*R2NH + CO2(g) <--> R2NCO2- + R2NH2+      

       REAC_DH(1)      = -71649.0    ! J/mol
       REAC_DS(1)      =  -200.0     ! J/mol.K
       REAC_E(1)       =  78728.0    ! J/mol
       REAC_LOGZETA(1) =  3.3115     ! [-] 

!   2. Wet CO2 Adsorption 
!      R2NH + H2O(g) + CO2(g) <--> HC03- + R2NH2+ 

       REAC_DH(2)      = -87733.0    ! J/mol
       REAC_DS(2)      =  -260.83    ! J/mol.K
       REAC_E(2)       =  11360.0    ! J/mol
       REAC_LOGZETA(2) =  0.6165     ! [-]

!   3. Water physical Adsorption 
!      H2O(g) <--> H2O(abs)   

       REAC_DH(3)      = -98828.0    ! J/mol
       REAC_DS(3)      =  -246.76    ! J/mol.K
       REAC_E(3)       =  67960.0    ! J/mol
       REAC_LOGZETA(3) =  4.2881     ! [-]

      RETURN  
      END SUBROUTINE USR_INIT_NAMELIST 
