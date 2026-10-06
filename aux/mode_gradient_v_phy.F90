MODULE MODE_GRADIENT_V_PHY

!$ACDC singlecolumn
!$ACDC manyblocks

IMPLICIT NONE
CONTAINS
     !     #######################################################
      SUBROUTINE GZ_V_VW_PHY(D,PA,PDZZ,PGZ_V_VW)
!     #######################################################
!
!!****  *GZ_V_VW - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the VW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     V point. The result is placed at the VW vorticity point.
!
!
!                   dzm(PA) 
!      PGZ_V_VW =   ------  
!                    ____y
!                    d*zz   
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MYM     : Shuman functions (mean operators)
!!      DZM     : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart        *INM and Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/07/94
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY, ONLY : DZM_PHY, MYM_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
!
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PGZ_V_VW ! result UW point
REAL, DIMENSION(D%NIJT,D%NKT) :: PA_WORK, PDZZ_WORK
!
INTEGER :: JIJ,JK
INTEGER :: IIJB,IIJE,IKT
!
!*       0.2   declaration of local variables
!
IIJE=D%NIJEC
IIJB=D%NIJB
IKT=D%NKT
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_V_VW
!              ---------------------
!
CALL DZM_PHY(D,PA,PA_WORK)
CALL MYM_PHY(D,PDZZ,PDZZ_WORK)
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PGZ_V_VW(JIJ, JK)= PA_WORK(JIJ, JK)/PDZZ_WORK(JIJ, JK)
  END DO
END DO

!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_V_VW_PHY
      SUBROUTINE GY_V_M_PHY(D,OFLAT,PA,PDYY,PDZZ,PDZY,PGY_V_M)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #######################################################
!
!!****  *GY_V_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     V point. The result is placed at the mass point.
!
!
!                       (          ______________z )
!                       (          (___________y ) )
!                    1  (          (d*zy dzm(PA) ) ) 
!      PGY_V_M =   ---- (dyf(PA) - (------------)) )
!                  ___y (          (             ) )
!                  d*yy (          (      d*zz   ) )     
!
!       
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MYF,MZF         : Shuman functions (mean operators)
!!      DYF,DZF         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart        *INM and Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    19/07/94
!!                  18/10/00 (V.Masson) add LFLAT switch
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY, ONLY : DZM_PHY, DYF_PHY, MYF_PHY, MZF_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
!
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PGY_V_M ! result mass point
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZWORK1, ZWORK2, ZWORK3, ZWORK4
INTEGER :: IIJB,IIJE,IKT
INTEGER :: JIJ,JK
!
!*       0.2   declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_V_M',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_V_M
!              --------------------
!
CALL DYF_PHY(D,PA,ZWORK1)
CALL MYF_PHY(D,PDYY,ZWORK2)
!
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PA,ZWORK3)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK3(JIJ, JK) = ZWORK3(JIJ, JK) * PDZY(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D,ZWORK3,ZWORK4)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK4(JIJ, JK) = ZWORK4(JIJ, JK) / PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MZF_PHY(D,ZWORK4,ZWORK3)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
        PGY_V_M(JIJ, JK) = ( ZWORK1(JIJ, JK) - ZWORK3(JIJ, JK)) &
                                           / ZWORK2(JIJ, JK)
    END DO
  END DO
  
ELSE
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGY_V_M(JIJ, JK)= ZWORK1(JIJ, JK) / ZWORK2(JIJ, JK)
    END DO
  END DO
  
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_V_M',1,ZHOOK_HANDLE)
END SUBROUTINE GY_V_M_PHY
!
END MODULE MODE_GRADIENT_V_PHY
