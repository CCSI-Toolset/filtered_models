!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments: Write reaction rates in units of moles/sec.cm^3 (cgs) or  !
!  kmoles/sec.m^3 (SI). Units should match those specified in the data !
!  file.
!                                                                      !
!  Example reaction: Methane combustion                                !
!                                                                      !
!  mfix.dat input:                                                     !
!``````````````````````````````````````````````````````````````````````!
!    @(RXNS)                                                           !
!      CH4_Comb { chem_eq = "CH4 + 2.0*O2 --> CO2 + 2.0*H2O" }         !
!    @(END)                                                            !
!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!  usr_rates.f input:                                                  !
!``````````````````````````````````````````````````````````````````````!
!    c_O2  = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))                          !
!    c_CH4 = (RO_g(IJK)*X_g(IJK,CH4)/MW_g(CH4))                        !
!    RATES(CH4_Comb) = 2.0d5 * EP_g(IJK) * c_O2 * c_CH4                !
!``````````````````````````````````````````````````````````````````````!
!  * Species alias and reaction names given in the data file can be    !
!    used in reference to the reaction index in RATES and a species    !
!    index in gas/solids phase variables.                              !
!                                                                      !
!  * Additional information is provided in section 5.11 of the code    !
!    Readme.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE sendrecv
      USE toleranc
      USE usr
      USE fun_avg

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
      real(kind=8), dimension(nReac) :: k,kappa,fwd,rev
      real(kind=8) :: C_CO2, C_H2O, p_CO2, p_H2O
      real(kind=8) :: C_R2NH, C_R2NCO2, C_R2NH2, C_H2O_s, C_HCO3,MW_MIX_gas
      real(kind=8) :: T_gas, T_solid

      integer :: iSpec, iReac, j, m , iSolid

! Reaction rates:
!`````````````````````````````````````````````````````````````````````//
! Include reaction rates here. Reaction rates should be stored in the
! variable RATES. The reaction name given in the data file can be used
! to store the rate in the appropriate array location. Additional
! input format parameters are given in Section 4.11 of the code Readme.

      RATES(:) = ZERO

!
! There are six reactions defined in mfix.dat, corresponding to three 
! forward and three reversed reactions, which parameters are defined in usr_mod:
!
!  1. iReacc = 1 : Dry CO2 Adsorption  2R2NH + CO2 <--> R2NH2 + R2NCO2
!                  IDs: 1 and 2
!
!                  ID  : 1
!                  Name: FWD_DRY_CO2_ADSORPTION
!                  Chemical Eq:  2R2NH + CO2 --> R2NH2 + R2NCO2
! 
!                  ID  : 2
!                  Name: REV_DRY_CO2_ADSORPTION
!                  Chemical Eq:  R2NH2 + R2NCO2 --> 2R2NH + CO2
!
!  2. iReacc = 2 : Wet CO2 Adsorption  R2NH + CO2 + H2O_S <--> R2NH2 + HCO3
!                  IDs: 3 and 4
!
!                  ID  : 3
!                  Name: FWD_WET_CO2_ADSORPTION
!                  Chemical Eq:  R2NH + CO2 + H2O_S --> R2NH2 + HCO3
!
!                  ID  : 4 
!                  Name: REV_WET_CO2_ADSORPTION
!                  Chemical Eq:  R2NH2 + HCO3 --> R2NH + CO2 + H2O_S
!
!  3. iReacc = 3 : H2O physisorption  H2O_G <--> H2O_S
!                  IDs: 5 and 6
!
!                  ID  : 5
!                  Name: FWD_H2O_PHYSISORPTION
!                  Chemical Eq:  H2O_G --> H2O_S
! 
!                  ID  : 6 
!                  Name: REV_H2O_PHYSISORPTION
!                  Chemical Eq:  H2O_S --> H2O_G



!
! Reactions occur betwee gas phase and solids phase 1
!
        iSolid = 1

!
!   Partial pressures
!
        p_CO2  = max(X_g(IJK,CO2),X_zero)   * P_G(IJK)  * (MW_MIX_g(IJK)/MW_g(CO2))  !modified by zhijie
        p_H2O  = max(X_g(IJK,H2O_g),X_zero) * P_G(IJK)  * (MW_MIX_g(IJK)/MW_g(H2O_g))  !modified by zhijie

!
!   Molar concentrations
!
        C_CO2  = RO_g(IJK) * max(X_g(IJK,CO2),X_zero) /MW_g(CO2)
        C_H2O  = RO_g(IJK) * max(X_g(IJK,H2O_g),X_zero) /MW_g(H2O_g)

        C_R2NH   = RO_s(IJK,iSolid) * max(X_s(IJK,iSolid,R2NH),X_zero)   / MW_s(iSolid,R2NH)
        C_R2NCO2 = RO_s(IJK,iSolid) * max(X_s(IJK,iSolid,R2NCO2),X_zero) / MW_s(iSolid,R2NCO2)
        C_R2NH2  = RO_s(IJK,iSolid) * max(X_s(IJK,iSolid,R2NH2),X_zero)  / MW_s(iSolid,R2NH2)
        C_HCO3   = RO_s(IJK,iSolid) * max(X_s(IJK,iSolid,HCO3),X_zero)   / MW_s(iSolid,HCO3)
        C_H2O_s  = RO_s(IJK,iSolid) * max(X_s(IJK,iSolid,H2O_s),X_zero)  / MW_s(iSolid,H2O_s)
!
!
!   Temperature of gas and solids phase (local IJK value)

        T_gas   = T_g(IJK)
        T_solid = T_s(IJK,iSolid)
!
!
!   k and kappa for each reaction
!
!
        if (use_conc) then
            do iReac = 1, nReac
                k(iReac) = reac(iReac)%A*(Rgas_press*T_gas)**reac(iReac)%m*T_solid*exp( -reac(iReac)%E/T_solid ) !Chris Montgomery - m exponent for RT
                kappa(iReac) = reac(iReac)%B*((Rgas_press*T_gas)/reac(iReac)%P_ref)*exp(reac(iReac)%dS - reac(iReac)%dH/T_solid )
            end do
        else
            do iReac = 1, nReac
                k(iReac) = reac(iReac)%A*T_solid*exp( -reac(iReac)%E/T_solid )
                kappa(iReac) = (reac(iReac)%B/reac(iReac)%P_ref)*exp(reac(iReac)%dS - reac(iReac)%dH/T_solid )
            end do
        end if

        k(3) = k(3)/1.e6   !  Artificially reduce k value for water physisorption



!
!   iReacc = 1 : Dry CO2 Adsorption
!
        iReac = 1

        if (use_conc) then
            fwd(iReac) = EP_S(IJK,iSolid)*k(iReac)*C_CO2**reac(iReac)%m*C_R2NH**2  !Chris Montgomery - m exponent for use_conc = .TRUE.
            rev(iReac) = EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_R2NCO2*C_R2NH2
        else
            fwd(iReac) = EP_S(IJK,iSolid)*k(iReac)* p_CO2**reac(iReac)%m * C_R2NH**2
            rev(iReac) = EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_R2NCO2*C_R2NH2
        end if

!     
!   iReacc = 2 : Wet CO2 Adsorption
!      
        iReac = 2
        if (use_conc) then
            fwd(iReac) =  EP_S(IJK,iSolid)*k(iReac)*C_CO2*C_R2NH*C_H2O_s
            rev(iReac) =  EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_R2NH2*C_HCO3
        else
            fwd(iReac) =  EP_S(IJK,iSolid)*k(iReac)*p_CO2*C_R2NH*C_H2O_s
            rev(iReac) =  EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_R2NH2*C_HCO3
        end if
      

!
!   iReacc = 3 : H2O physisorption
!
        iReac = 3
        if (use_conc) then
             fwd(iReac) = EP_S(IJK,iSolid)*k(iReac)*C_H2O
             rev(iReac) = EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_H2O_s
        else
             fwd(iReac) = EP_S(IJK,iSolid)*k(iReac)*p_H2O
             rev(iReac) = EP_S(IJK,iSolid)*k(iReac)/kappa(iReac)*C_H2O_s
        end if

!
!   Store reaction rates
!

      RATES(fwd_Dry_CO2_Adsorption)  = fwd(1)
      RATES(rev_Dry_CO2_Adsorption)  = rev(1)
      RATES(fwd_Wet_CO2_Adsorption)  = fwd(2)
      RATES(rev_Wet_CO2_Adsorption)  = rev(2)
      RATES(fwd_H2O_physisorption)   = fwd(3)
      RATES(rev_H2O_physisorption)   = rev(3)



      if (nRR .ge. 1) ReactionRates(ijk, fwd_Dry_CO2_Adsorption) = RATES(fwd_Dry_CO2_Adsorption)
      if (nRR .ge. 2) ReactionRates(ijk, rev_Dry_CO2_Adsorption) = RATES(rev_Dry_CO2_Adsorption)
      if (nRR .ge. 3) ReactionRates(ijk, fwd_Wet_CO2_Adsorption) = RATES(fwd_Wet_CO2_Adsorption)
      if (nRR .ge. 4) ReactionRates(ijk, rev_Wet_CO2_Adsorption) = RATES(rev_Wet_CO2_Adsorption)
      if (nRR .ge. 5) ReactionRates(ijk, fwd_H2O_physisorption)  = RATES(fwd_H2O_physisorption) 
      if (nRR .ge. 6) ReactionRates(ijk, rev_H2O_physisorption)  = RATES(rev_H2O_physisorption) 



      RETURN  

      END SUBROUTINE USR_RATES
