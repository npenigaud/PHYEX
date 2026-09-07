!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    #######################  
     MODULE MODE_ROTATE_WIND

!$ACDC singlecolumn --dummy
!$ACDC manyblocks --dummy

!    #######################
IMPLICIT NONE
CONTAINS
!     ###########################################################
      SUBROUTINE ROTATE_WIND(D,PU,PV,PW,                          &
                             PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,   &
                             PCOSSLOPE,PSINSLOPE,               &
                             PDXX,PDYY,PDZZ,                    &
                             PUSLOPE,PVSLOPE                    )
!     ###########################################################
!
!
!!****  *ROTATE_WIND* - computes the wind components along the maximum slope 
!!               direction and its normal direction in the first mass level.
!!
!!    PURPOSE
!!    -------
!!**** 
!        The purpose of this routine is to compute the wind component parallel 
!     to the orography at the first mass level. The exact location where these
!     components are computed is the point of intersection between the normal 
!     to the orography and the first mass-level hyper-plane at PDZZ(:,:,IKB)/2 
!        
!!**  METHOD
!!    ------
!!       The values of the 3 cartesian components of the wind are determined
!!    by a bilinear interpolation between the 4 nearest points in the first 
!!    mass-level hyper-plane. These points are found according to the signs of 
!!    the slopes' sinus and cosinus. For each direction of interpolation, the 
!!    two different localizations (mass or flux grids) are used to avoid 
!!    lateral boundary problems.  
!!       Then, the rotation is performed for the wind components. The rotation 
!!    angle is the angle between the x axe and the maximum slope direction 
!!    defined by the slope vector (dZs/dx , dZs/dy).
!!        Finally, the horizontal components are set at the marginal points 
!!    according to cyclic boundary conditions because this is the only case
!!    where these points can be considered.
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS 
!!    ------------------
!!
!!       MODD_CONF      : L2D   switch for 2D model version
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joel Stein              * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         14/11/95
!!      Modifications:   15/05/96, (N. wood)
!!                                 take into account no slip conditions 
!!                                 at the surface
!!                       14/02/01  (V. Masson)
!!                                 Slip condition at the surface restored
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
! These macro are handled by pft_tool.py --craybyPassDOCONCURRENT applied on Cray Rules
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN) :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PU,PV,PW        ! cartesian components
                                 ! of the wind
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDXX, PDYY, PDZZ
                                 ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(OUT)  ::  PUSLOPE         ! wind component along 
                                 ! the maximum slope direction
REAL, DIMENSION(D%NIJT),   INTENT(OUT)  ::  PVSLOPE         ! wind component along
                                 !  the direction normal to the maximum slope one
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(D%NIJT) :: ILOC,JLOC
              ! shift index to find the 4 nearest points in x and y directions
REAL,    DIMENSION(D%NIJT) :: ZCOEFF,ZCOEFM,     &
              ! interpolation weigths for flux and mass locations
                                                           ZUINT,ZVINT,ZWINT, &
              ! intermediate values of the cartesian components after x interp.
                                                           ZUFIN,ZVFIN,ZWFIN, &
              ! final values of the cartesian components after the 2 interp.
                                                           ZWGROUND
              ! vertical velocity at the surface                                                            
INTEGER     :: IIB,IIE,IJB,IJE,IKB
              ! index values for the Beginning or the End of the physical 
              ! domain in x,y and z directions
INTEGER     :: IIU,IJU
              ! arrays' sizes for i and j indices
INTEGER     :: JI,JJ
!      
!----------------------------------------------------------------------------
!
!*      1.    PRELIMINARIES
!             -------------
!
!
IIB = 2
IJB = 2
IIU = SIZE(PU,1)
IJU = SIZE(PU,2)
IIE = IIU - 1
IJE = IJU - 1
IKB = D%NKB
!



!!provisoire pour changement interface!!!!DO JJ=1, IJU
!!provisoire pour changement interface!!!!  DO JI=1, IIU
!!provisoire pour changement interface!!!!    PUSLOPE(JI, JJ)=0.
!!provisoire pour changement interface!!!!    PVSLOPE(JI, JJ)=0.
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    ZWGROUND(JI, JJ) = PW(JI, JJ, IKB)
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    !*      2.    INTERPOLATE THE CARTESIAN COMPONENTS
!!provisoire pour changement interface!!!!    !             ------------------------------------
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    ILOC(JI, JJ)=NINT(SIGN(1.,-PCOSSLOPE(JI, JJ)))
!!provisoire pour changement interface!!!!    JLOC(JI, JJ)=NINT(SIGN(1.,-PSINSLOPE(JI, JJ)))
!!provisoire pour changement interface!!!!  END DO
!!provisoire pour changement interface!!!!END DO
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!! interpolation in x direction
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!DO JJ=1,IJU
!!provisoire pour changement interface!!!!  DO JI=IIB,IIE
!!provisoire pour changement interface!!!!    ZCOEFF(JI,JJ) =                                                  &
!!provisoire pour changement interface!!!!      (0.5*PDXX(JI,JJ,IKB) + 0.5*PDZZ(JI,JJ,IKB)*PDIRCOSXW(JI,JJ) )  & 
!!provisoire pour changement interface!!!!      * 2. / (PDXX(JI,JJ,IKB)+PDXX(JI+1,JJ,IKB))
!!provisoire pour changement interface!!!!    ZUINT(JI,JJ) = ZCOEFF(JI,JJ)      * PU(JI+1,JJ,IKB)  +           &
!!provisoire pour changement interface!!!!                   (1.-ZCOEFF(JI,JJ)) * PU(JI,JJ,IKB)
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSXW(JI,JJ))      & 
!!provisoire pour changement interface!!!!                             / PDXX(JI+(ILOC(JI,JJ)+1)/2,JJ,IKB)
!!provisoire pour changement interface!!!!    ZVINT(JI,JJ) = ZCOEFM(JI,JJ)      * PV(JI,JJ,IKB)              +        &
!!provisoire pour changement interface!!!!                   (1.-ZCOEFM(JI,JJ)) * PV(JI+ILOC(JI,JJ),JJ,IKB)
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    ZWINT(JI,JJ) = ZCOEFM(JI,JJ)  * (PW(JI,JJ,IKB+1)+ZWGROUND(JI,JJ)) * 0.5    &
!!provisoire pour changement interface!!!!              + (1.-ZCOEFM(JI,JJ))                                             &
!!provisoire pour changement interface!!!!               *(PW(JI+ILOC(JI,JJ),JJ,IKB+1)+ZWGROUND(JI+ILOC(JI,JJ),JJ)) * 0.5
!!provisoire pour changement interface!!!!  END DO
!!provisoire pour changement interface!!!!END DO
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!! interpolation in y direction
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!DO JJ=IJB,IJE
!!provisoire pour changement interface!!!!  DO JI=IIB,IIE
!!provisoire pour changement interface!!!!    ZCOEFF(JI,JJ) =                                                     &
!!provisoire pour changement interface!!!!      (0.5*PDYY(JI,JJ,IKB) + 0.5*PDZZ(JI,JJ,IKB)*PDIRCOSYW(JI,JJ) )     & 
!!provisoire pour changement interface!!!!      * 2. / (PDYY(JI,JJ,IKB)+PDYY(JI+1,JJ,IKB))
!!provisoire pour changement interface!!!!    ZVFIN(JI,JJ) = ZCOEFF(JI,JJ)      * ZVINT(JI,JJ+1)  +               &
!!provisoire pour changement interface!!!!                   (1.-ZCOEFF(JI,JJ)) * ZVINT(JI,JJ)
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSYW(JI,JJ))   & 
!!provisoire pour changement interface!!!!                             / PDYY(JI,JJ+(JLOC(JI,JJ)+1)/2,IKB)
!!provisoire pour changement interface!!!!    ZUFIN(JI,JJ) = ZCOEFM(JI,JJ)      * ZUINT(JI,JJ)                +    &
!!provisoire pour changement interface!!!!                   (1.-ZCOEFM(JI,JJ)) * ZUINT(JI,JJ+JLOC(JI,JJ))
!!provisoire pour changement interface!!!!    ZWFIN(JI,JJ) = ZCOEFM(JI,JJ)      * ZWINT(JI,JJ)                +    &
!!provisoire pour changement interface!!!!                   (1.-ZCOEFM(JI,JJ)) * ZWINT(JI,JJ+JLOC(JI,JJ))
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!!*      3.    ROTATE THE WIND
!!provisoire pour changement interface!!!!!             ---------------
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!    PUSLOPE(JI,JJ) = PCOSSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZUFIN(JI,JJ) +   &
!!provisoire pour changement interface!!!!                     PSINSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZVFIN(JI,JJ) +   &
!!provisoire pour changement interface!!!!                            SQRT(1.-PDIRCOSZW(JI,JJ)**2) * ZWFIN(JI,JJ)
!!provisoire pour changement interface!!!!    !              
!!provisoire pour changement interface!!!!    PVSLOPE(JI,JJ) =-PSINSLOPE(JI,JJ)                    * ZUFIN(JI,JJ) +   &
!!provisoire pour changement interface!!!!                     PCOSSLOPE(JI,JJ)                    * ZVFIN(JI,JJ)
!!provisoire pour changement interface!!!!    !
!!provisoire pour changement interface!!!!  END DO
!!provisoire pour changement interface!!!!END DO
!

!
!
!----------------------------------------------------------------------------
!
END SUBROUTINE ROTATE_WIND
!
!     ##############################################
      SUBROUTINE UPDATE_ROTATE_WIND(D,PUSLOPE,PVSLOPE,HLBCX,HLBCY)
!     ##############################################
!!
!!****  *UPDATE_ROTATE_WIND* routine to set rotate wind values at the border
!
!!    AUTHOR
!!    ------
!!
!!     P Jabouille   *CNRM METEO-FRANCE
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   24/06/99
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll, ONLY: ADD2DFIELD_ll, UPDATE_HALO_ll, CLEANLIST_ll, &
                   LWEST_ll, LEAST_ll, LSOUTH_ll, LNORTH_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
#ifdef MNH_OPENACC
USE MODI_GET_HALO, ONLY: GET_2D_HALO_DDC
#endif
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),       INTENT(IN) :: D
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PUSLOPE,PVSLOPE
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
logical                :: gwest, geast, gnorth, gsouth
INTEGER :: IIB, IIE, IJB, IJE
!
! tangential surface fluxes in the axes following the orography
!
!*        1  PROLOGUE
!
NULLIFY(TZFIELDS_ll)
!
IIB = D%NIB
IJB = D%NJB
IIE = D%NIE
IJE = D%NJE
!         2 Update halo if necessary
!
!!provisoire pour changement interface!!!!!!!$IF (NHALO == 1) THEN
!!provisoire pour changement interface!!!!!#ifndef MNH_OPENACC
!!provisoire pour changement interface!!!!!  CALL ADD2DFIELD_ll( TZFIELDS_ll, PUSLOPE, 'UPDATE_ROTATE_WIND::PUSLOPE' )
!!provisoire pour changement interface!!!!!  CALL ADD2DFIELD_ll( TZFIELDS_ll, PVSLOPE, 'UPDATE_ROTATE_WIND::PVSLOPE' )
!!provisoire pour changement interface!!!!!  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
!!provisoire pour changement interface!!!!!  CALL CLEANLIST_ll(TZFIELDS_ll)
!!provisoire pour changement interface!!!!!#else
!!provisoire pour changement interface!!!!!!
!!provisoire pour changement interface!!!!!!  /!\ warning conner needed -> GET_HALO...C
!!provisoire pour changement interface!!!!!!
!!provisoire pour changement interface!!!!!CALL GET_2D_HALO_DDC( PUSLOPE, HNAME='UPDATE_ROTATE_WIND::PUSLOPE' )
!!provisoire pour changement interface!!!!!CALL GET_2D_HALO_DDC( PVSLOPE, HNAME='UPDATE_ROTATE_WIND::PVSLOPE' )
!!provisoire pour changement interface!!!!!#endif
!!provisoire pour changement interface!!!!!!!$ENDIF
!!provisoire pour changement interface!!!!!!
!!provisoire pour changement interface!!!!!!        3 Boundary conditions for non cyclic case
!!provisoire pour changement interface!!!!!!
!!provisoire pour changement interface!!!!!gwest  = HLBCX(1) /= "CYCL" .AND. LWEST_ll()
!!provisoire pour changement interface!!!!!geast  = HLBCX(2) /= "CYCL" .AND. LEAST_ll()
!!provisoire pour changement interface!!!!!gsouth = HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()
!!provisoire pour changement interface!!!!!gnorth = HLBCY(2) /= "CYCL" .AND. LNORTH_ll()
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!!
!!provisoire pour changement interface!!!!!IF ( gwest ) THEN
!!provisoire pour changement interface!!!!!  PUSLOPE(IIB-1,:)=PUSLOPE(IIB,:)
!!provisoire pour changement interface!!!!!  PVSLOPE(IIB-1,:)=PVSLOPE(IIB,:)
!!provisoire pour changement interface!!!!!END IF
!!provisoire pour changement interface!!!!!IF ( geast ) THEN
!!provisoire pour changement interface!!!!!  PUSLOPE(IIE+1,:)=PUSLOPE(IIE,:)
!!provisoire pour changement interface!!!!!  PVSLOPE(IIE+1,:)=PVSLOPE(IIE,:)
!!provisoire pour changement interface!!!!!END IF
!!provisoire pour changement interface!!!!!IF ( gsouth ) THEN
!!provisoire pour changement interface!!!!!  PUSLOPE(:,IJB-1)=PUSLOPE(:,IJB)
!!provisoire pour changement interface!!!!!  PVSLOPE(:,IJB-1)=PVSLOPE(:,IJB)
!!provisoire pour changement interface!!!!!END IF
!!provisoire pour changement interface!!!!!IF( gnorth ) THEN
!!provisoire pour changement interface!!!!!  PUSLOPE(:,IJE+1)=PUSLOPE(:,IJE)
!!provisoire pour changement interface!!!!!  PVSLOPE(:,IJE+1)=PVSLOPE(:,IJE)
!!provisoire pour changement interface!!!!!END IF
!
END SUBROUTINE UPDATE_ROTATE_WIND
END MODULE MODE_ROTATE_WIND
