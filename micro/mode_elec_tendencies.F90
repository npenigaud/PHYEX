!MNH_LIC Copyright 2022-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ELEC_TENDENCIES
!
IMPLICIT NONE
CONTAINS
!
!     #########################################################################################
      SUBROUTINE ELEC_TENDENCIES (D, CST, ICED, ICEP, ELECD, ELECP,                           &
                                  KRR, KMICRO, PTSTEP, ODMICRO,                               &
                                  BUCONF, TBUDGETS, KBUDGETS,                                 &
                                  HCLOUD, PTHVREFZIKB,                                        &
                                  PRHODREF, PRHODJ, PZT, PCIT,                                &
                                  PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                         &
                                  PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                 &
                                  PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                 &
                                  PRVHENI, PRRHONG, PRIMLTC,                                  &
                                  PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG,                &
                                  PRCAUTR, PRCACCR, PRREVAV,                                  &
                                  PRCRIMSS, PRCRIMSG, PRSRIMCG, PRRACCSS, PRRACCSG, PRSACCRG, &
                                  PRSMLTG, PRICFRRG, PRRCFRIG,                                &
                                  PRCWETG, PRIWETG, PRRWETG, PRSWETG,                         &
                                  PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG,                         &
                                  PRGMLTR, PRCBERI,                                           & 
                                  PRCMLTSR, PRICFRR,                                          & !- opt. param. for ICE3 
                                  PCCT, PCRT, PCST, PCGT,                                     & !-- optional
                                  PRVHENC, PRCHINC, PRVHONH,                                  & !| parameters 
                                  PRRCVRC, PRICNVI, PRVDEPI, PRSHMSI, PRGHMGI,                & !|    for
                                  PRICIBU, PRIRDSF,                                           & !|    LIMA
                                  PRCCORR2, PRRCORR2, PRICORR2,                               & !--
                                  PRWETGH, PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH,       & !--  optional
                                  PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH,                & !|  parameters
                                  PRHMLTR, PRDRYHG,                                           & !|     for
                                  PRHT, PCHT, PQHT, PQHS)                                 !--    hail

!$ACDC singlecolumn --dummy
!$ACDC manyblocks --dummy --ignore-array-variables=PBUDGETS

!     ##########################################################################################
!
!!****  * - compute the explicit cloud electrification sources
!!
!!    This routine is adapted from rain_ice_elec.f90.
!!    To avoid duplicated routines, the cloud electrification routine is now CALLed 
!!    at the end of the microphysics scheme but needs the microphysical tendencies as arguments.
!!    The sedimentation source for electric charges is treated separately.
!!
!!    AUTHOR
!!    ------
!!      C. Barthe    * LAERO *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    February 2022
!!
!!      Modifications
!! C. Barthe   12/04/2022   include electrification from LIMA
!! C. Barthe   22/03/2023   5-6: take into account news from LIMA (Ns, Ng, Nh, CIBU and RDSF) and PHYEX
!! C. Barthe   13/07/2023   5-6: Ns, Ng and Nh can be pronostic variables (LIMA2)
!! C. Barthe   22/11/2023   initialize Nx to 0 when 1-moment
!!
!------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_BUDGET,     ONLY:  NBUDGET_SV1, TBUDGETDATA_PTR, TBUDGETCONF_T
!
USE MODD_CST,              ONLY: CST_t
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_NSV, ONLY: NSV_ELECBEG, NSV_ELECEND
USE MODD_ELEC_DESCR, ONLY: CNI_CHARGING, ELEC_DESCR_T, LINDUCTIVE, XFG, XFR, XFS, XLBDAG_MAXE, &
                           XLBDAS_MAXE, XLIM_NI_IG, XLIM_NI_IS, XLIM_NI_SG, XQTC, XQTMIN, XRTMIN_ELEC      
USE MODD_ELEC_PARAM, ONLY: ELEC_PARAM_T, NIND_LWC, NIND_TEMP, XAIGAMMABI, XAUX_LIM, XAUX_LIM1, &
                           XAUX_LIM2, XAUX_LIM3, XCOEF_RQ_C, XCOEF_RQ_G, XCOEF_RQ_H, XCOEF_RQ_I, &
                           XCOEF_RQ_R, XCOEF_RQ_S, XESMIN, XEXQRCFRIG, XEXQSRIMCG, XFQIAGGSBG, &
                           XFQIAGGSBH, XFQIAGGSBS, XFQIAGGSBT1, XFQIAGGSBT2, XFQIAGGSBT3, XFQIAGGSN, &
                           XFQIAGGSN_TAK, XFQIAGGSP, XFQIAGGSP_TAK, XFQIDRYGBG, XFQIDRYGBS, &
                           XFQIDRYGBSN, XFQIDRYGBSN_TAK, XFQIDRYGBSP, XFQIDRYGBSP_TAK, XFQIDRYGBT1, &
                           XFQIDRYGBT2, XFQIDRYGBT3, XFQRACCS, XFQRDRYG, XFQSDRYG, XFQSDRYGBG, &
                           XFQSDRYGBH, XFQSDRYGBS, XFQSDRYGBT1, XFQSDRYGBT10, XFQSDRYGBT11, &
                           XFQSDRYGBT2, XFQSDRYGBT3, XFQSDRYGBT4, XFQSDRYGBT5, XFQSDRYGBT6, &
                           XFQSDRYGBT7, XFQSDRYGBT8, XFQSDRYGBT9, XGAMINC_RIM3, XHIDRYG, XIMN, &
                           XIMP, XIND1, XIND2, XIND3, XINN, XINP, XKER_Q_LIMSG, XKER_Q_RACCS, &
                           XKER_Q_RACCSS, XKER_Q_RDRYG, XKER_Q_SACCRG, XKER_Q_SDRYG, XKER_Q_SDRYGB, &
                           XKER_Q_SDRYGB1, XKER_Q_SDRYGB2, XLBQIDRYGBG, XLBQRACCS1, XLBQRACCS2, &
                           XLBQRACCS3, XLBQRDRYG1, XLBQRDRYG2, XLBQRDRYG3, XLBQSACCRG1, XLBQSACCRG2, &
                           XLBQSACCRG3, XLBQSDRYG1, XLBQSDRYG2, XLBQSDRYG3, XLBQSDRYGB1SN, &
                           XLBQSDRYGB1SP, XLBQSDRYGB2SN, XLBQSDRYGB2SP, XLBQSDRYGB3SN, XLBQSDRYGB3SP, &
                           XLBQSDRYGB4G, XLBQSDRYGB4H, XLBQSDRYGB5G, XLBQSDRYGB5H, XLBQSDRYGB6G, &
                           XLBQSDRYGB6H, XLWCC, XMANSELL, XQHON, XQRCFRIG, XQSRIMCG, XSAUNDER, &
                           XSKN, XSKN_TAK, XSKP, XSKP_TAK, XSMN, XSMP, XSNN, XSNP, XTAKA_TM, &
                           XVGCOEF, XVSCOEF      
USE MODD_ELEC_n, ONLY: XEFIELDW, XEW, XIND_RATE, XNI_IAGGS, XNI_IDRYG, XNI_SDRYG
USE MODD_PARAM_LIMA,       ONLY: XALPHAI_L=>XALPHAI, XNUI_L=>XNUI,   &
                                 XCEXVT_L=>XCEXVT, XRTMIN_L=>XRTMIN, &
                                 LCIBU, LRDSF,                       &
                                 NMOM_C, NMOM_R, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD,  ONLY: XAI_L=>XAI, XBI_L=>XBI,   &
                                 XDS_L=>XDS, XCXS_L=>XCXS, &
                                 XCOLEXIS_L=>XCOLEXIS
USE MODD_PARAM_LIMA_MIXED, ONLY: XDG_L=>XDG, XCXG_L=>XCXG,                           &
                                 XCOLIG_L=>XCOLIG, XCOLEXIG_L=>XCOLEXIG,             &
                                 XCOLSG_L=>XCOLSG, XCOLEXSG_L=>XCOLEXSG,             &
                                 NGAMINC_L=>NGAMINC,                                 &
                                 NACCLBDAR_L=>NACCLBDAR, NACCLBDAS_L=>NACCLBDAS,     &
                                 XACCINTP1S_L=>XACCINTP1S, XACCINTP2S_L=>XACCINTP2S, &
                                 XACCINTP1R_L=>XACCINTP1R, XACCINTP2R_L=>XACCINTP2R, &
                                 NDRYLBDAR_L=>NDRYLBDAR, NDRYLBDAS_L=>NDRYLBDAS,     &
                                 NDRYLBDAG_L=>NDRYLBDAG,                             &
                                 XDRYINTP1R_L=>XDRYINTP1R, XDRYINTP2R_L=>XDRYINTP2R, &
                                 XDRYINTP1S_L=>XDRYINTP1S, XDRYINTP2S_L=>XDRYINTP2S, &
                                 XDRYINTP1G_L=>XDRYINTP1G, XDRYINTP2G_L=>XDRYINTP2G, &
                                 XRIMINTP1_L=>XRIMINTP1,   XRIMINTP2_L=>XRIMINTP2
!
!#ifdef MNH_PGI
!USE MODE_PACK_PGI
!#endif
use mode_tools,           only: Countjv
!
USE MODE_COMPUTE_LAMBDA, ONLY: COMPUTE_LAMBDA
USE MODE_ELEC_COMPUTE_EX,ONLY: ELEC_COMPUTE_EX
USE MODI_MOMG, ONLY: MOMG
!
IMPLICIT NONE
!
!
!*      0.1   Declaration of dummy arguments
!
TYPE(DIMPHYEX_t),                      INTENT(IN)     :: D
TYPE(CST_t),                           INTENT(IN)     :: CST
TYPE(TBUDGETCONF_t),                   INTENT(IN)     :: BUCONF        ! budget structure
TYPE(RAIN_ICE_DESCR_t),                INTENT(IN)     :: ICED
TYPE(RAIN_ICE_PARAM_t),                INTENT(IN)     :: ICEP
TYPE(ELEC_PARAM_t),                    INTENT(IN)     :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),                    INTENT(IN)     :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS),INTENT(INOUT)  :: TBUDGETS
INTEGER,                               INTENT(IN)     :: KBUDGETS
!
INTEGER,                              INTENT(IN)    :: KMICRO
REAL,                                 INTENT(IN)    :: PTSTEP  ! Double Time step
                                                               ! (single if cold start)
INTEGER,                              INTENT(IN)    :: KRR     ! Number of moist variable
!
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PZT     ! Temperature (K)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCIT    ! Pristine ice n.c. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQPIT   ! Positive ion (Nb/kg) at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQNIT   ! Negative ion (Nb/kg) at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQCT    ! Cloud water charge at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQRT    ! Raindrops charge at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQIT    ! Pristine ice charge at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQST    ! Snow/aggregates charge at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PQGT    ! Graupel charge at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQPIS   ! Positive ion (Nb/kg) source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQNIS   ! Negative ion (Nb/kg) source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQCS    ! Cloud water charge source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQRS    ! Raindrops charge source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQIS    ! Pristine ice charge source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQSS    ! Snow/aggregates charge source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PQGS    ! Graupel charge source
!
! microphysics rates common to ICE3 and LIMA
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVHENI, &  ! heterogeneous nucleation mixing ratio change (HIND for LIMA)
                                                       PRCHONI, &  ! Homogeneous nucleation
                                                       PRRHONG, &  ! Spontaneous freezing mixing ratio change
                                                       PRVDEPS, &  ! Deposition on r_s,
                                                       PRIAGGS, &  ! Aggregation on r_s
                                                       PRIAUTS, &  ! Autoconversion of r_i for r_s production (CNVS for LIMA)
                                                       PRVDEPG, &  ! Deposition on r_g
                                                       PRCAUTR, &  ! Autoconversion of r_c for r_r production
                                                       PRCACCR, &  ! Accretion of r_c for r_r production
                                                       PRREVAV, &  ! Evaporation of r_r
                                                       PRIMLTC, &  ! Cloud ice melting mixing ratio change
                                                       PRCBERI, &  ! Bergeron-Findeisen effect
                                                       PRSMLTG, &  ! Conversion-Melting of the aggregates
                                                       PRRACCSS, PRRACCSG, PRSACCRG, & ! Rain accretion onto the aggregates
                                                       PRCRIMSS, PRCRIMSG, PRSRIMCG, & ! Cloud droplet riming of the aggregates
                                                       PRICFRRG, PRRCFRIG,           & ! Rain contact freezing
                                                       PRCWETG, PRIWETG, PRRWETG, PRSWETG, &  ! Graupel wet growth
                                                       PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, &  ! Graupel dry growth
                                                       PRGMLTR     ! Melting of the graupel
! microphysics rates specific to ICE3 (knmoments==1)
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRCMLTSR,&  ! Cld droplet collection onto aggregates by pos. temp.
                                                               PRICFRR     ! Rain contact freezing (part of ice crystals converted to rain)
! microphysics rates specific to LIMA (knmoments==2)
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRVHENC, &  ! Cld droplet formation
                                                               PRCHINC, &  ! Heterogeneous nucleation of coated IFN
                                                               PRVHONH, &  ! Nucleation of haze
                                                               PRRCVRC, &  ! Conversion of small drops into droplets
                                                               PRICNVI, &  ! Conversion snow --> ice
                                                               PRVDEPI, &  ! Deposition on r_i
                                                               PRSHMSI, PRGHMGI, & ! Hallett Mossop for snow and graupel
                                                               PRICIBU, &  ! Collisional ice breakup
                                                               PRIRDSF, &  ! Raindrop shattering by freezing
                                                               PRCCORR2, PRRCORR2, PRICORR2 ! Correction inside LIMA splitting
! microphysics rates related to hail (krr == 7, lhail = .t.)
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRWETGH, &  ! Conversion of graupel into hail
                                                               PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, & ! Dry growth of hail
                                                               PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, & ! Wet growth of hail
                                                               PRHMLTR, &                                     ! Melting of hail  
                                                               PRDRYHG     ! Conversion of hail into graupel
!
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCCT   ! Cloud droplets conc. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCRT   ! Raindrops conc. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCST   ! Snow conc. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCGT   ! Graupel conc. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCHT   ! Hail conc. at t
!
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRHT   ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PQHT   ! Hail charge at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(INOUT) :: PQHS   ! Hail charge source
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!
!
!*      0.2   Declaration of local variables    
!
INTEGER :: II, JL         ! Loop indexes
INTEGER :: IIB, IIE, &    ! Define the domain
           IJB, IJE, &    ! where the microphysical sources
           IKB, IKE       ! must be computed
INTEGER                    :: IMICRO   ! nb of pts where r_x > 0
INTEGER, DIMENSION(KMICRO) :: I1
INTEGER, DIMENSION(KMICRO) :: II1, II2, II3
!
LOGICAL, DIMENSION(KMICRO) :: GMASK    !          Mask 
!REAL,    DIMENSION(KMICRO) :: ZMASK    !       to reduce
INTEGER                    :: IGMASK   ! the computation domain
!
REAL, DIMENSION(KMICRO) :: ZRHODREF  ! Reference density
REAL, DIMENSION(KMICRO) :: ZRHODJ    ! RHO times Jacobian
REAL, DIMENSION(KMICRO) :: ZZT     ! Temperature
!
REAL, DIMENSION(KMICRO) :: ZRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KMICRO) :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KMICRO) :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(KMICRO) :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(KMICRO) :: ZCCT    ! Cloud water conc. at t
REAL, DIMENSION(KMICRO) :: ZCRT    ! Raindrops conc. at t
REAL, DIMENSION(KMICRO) :: ZCIT    ! Pristine ice conc. at t
REAL, DIMENSION(KMICRO) :: ZCST    ! Snow/aggregate conc. at t
REAL, DIMENSION(KMICRO) :: ZCGT    ! Graupel conc. at t
REAL, DIMENSION(KMICRO) :: ZCHT    ! Hail conc. at t
!
REAL, DIMENSION(KMICRO) :: ZQPIT   ! Positive ion (/kg) at t
REAL, DIMENSION(KMICRO) :: ZQNIT   ! Negative ion (/kg) at t
REAL, DIMENSION(KMICRO) :: ZQCT    ! Cloud water charge at t
REAL, DIMENSION(KMICRO) :: ZQRT    ! Raindrops charge at t
REAL, DIMENSION(KMICRO) :: ZQIT    ! Pristine ice charge at t
REAL, DIMENSION(KMICRO) :: ZQST    ! Snow/aggregate charge at t
REAL, DIMENSION(KMICRO) :: ZQGT    ! Graupel charge at t
REAL, DIMENSION(KMICRO) :: ZQHT    ! Hail charge at t
!
REAL, DIMENSION(KMICRO) :: ZQPIS   ! Positive ion (/kg) source
REAL, DIMENSION(KMICRO) :: ZQNIS   ! Negative ion (/kg) source
REAL, DIMENSION(KMICRO) :: ZQCS    ! Cloud water charge source
REAL, DIMENSION(KMICRO) :: ZQRS    ! Raindrops charge source
REAL, DIMENSION(KMICRO) :: ZQIS    ! Pristine ice charge source
REAL, DIMENSION(KMICRO) :: ZQSS    ! Snow/aggregate charge source
REAL, DIMENSION(KMICRO) :: ZQGS    ! Graupel charge source
REAL, DIMENSION(KMICRO) :: ZQHS    ! Hail charge source
!
REAL, DIMENSION(KMICRO) :: ZLBDAC  ! Slope parameter of the droplets distribution
REAL, DIMENSION(KMICRO) :: ZLBDAR  ! Slope parameter of the raindrop distribution
REAL, DIMENSION(KMICRO) :: ZLBDAI  ! Slope parameter of the pristine ice distribution
REAL, DIMENSION(KMICRO) :: ZLBDAS  ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KMICRO) :: ZLBDAG  ! Slope parameter of the graupel distribution
REAL, DIMENSION(KMICRO) :: ZLBDAH  ! Slope parameter of the hail distribution
!
REAL, DIMENSION(KMICRO) :: ZECT    !
REAL, DIMENSION(KMICRO) :: ZERT    !     e_x coef
REAL, DIMENSION(KMICRO) :: ZEIT    !      in the 
REAL, DIMENSION(KMICRO) :: ZEST    ! q_x - D_x relation
REAL, DIMENSION(KMICRO) :: ZEGT    !
REAL, DIMENSION(KMICRO) :: ZEHT    !
!
LOGICAL, DIMENSION(KMICRO,4) :: GELEC  ! Mask for non-inductive charging
!
REAL, DIMENSION(:), ALLOCATABLE :: ZDQ, ZDQ_IS, ZDQ_IG, ZDQ_SG
!
! Non-inductive charging process following Gardiner et al. (1995)
REAL, DIMENSION(:), ALLOCATABLE :: ZDELTALWC ! Gap between LWC and a critical LWC
REAL, DIMENSION(:), ALLOCATABLE :: ZFT       ! Fct depending on temperature
!
! Non-inductive charging process following Saunders et al. (1991) / EW
REAL, DIMENSION(:), ALLOCATABLE :: ZEW      ! Effective liquid water content
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK  ! constant B 
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM  !  d_s exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN  !  v_g-v_s
REAL, DIMENSION(:), ALLOCATABLE :: ZFQIAGGS, ZFQIDRYGBS
REAL, DIMENSION(:), ALLOCATABLE :: ZLBQSDRYGB1S, ZLBQSDRYGB2S, ZLBQSDRYGB3S
!
! Non-inductive charging process following Saunders and Peck (1998) / RAR
REAL, DIMENSION(:), ALLOCATABLE :: ZRAR        ! Rime accretion rate
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IS  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IS  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IG  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IG  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK_SG  ! constant B 
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM_SG  !  d_s exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN_SG  !  v_g-v_s
!
! Inductive charging process (Ziegler et al., 1991)
REAL, DIMENSION(:), ALLOCATABLE :: ZEFIELDW  ! Vertical component of the electric field
!
REAL, DIMENSION(KMICRO) :: ZLIMIT   ! Used to limit the charge separated during NI process 
REAL, DIMENSION(KMICRO) :: ZQCOLIS  ! Collection efficiency between ice and snow
REAL, DIMENSION(KMICRO) :: ZQCOLIG  ! Collection efficiency between ice and graupeln
REAL, DIMENSION(KMICRO) :: ZQCOLSG  ! Collection efficiency between snow and graupeln
!
REAL                    :: ZRHO00, ZCOR00   ! Surface reference air density
REAL, DIMENSION(KMICRO) :: ZRHOCOR  ! Density correction for fallspeed
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1, IVEC2                   ! Vectors of indices for interpolation
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC2                   ! Work vectors for interpolation
REAL,    DIMENSION(:), ALLOCATABLE :: ZVECQ1, ZVECQ2, ZVECQ3, ZVECQ4 ! Work vectors for interpolation
!
REAL,    DIMENSION(KMICRO)   :: ZWQ, ZWQ_NI                    !  Work arrays
REAL,    DIMENSION(KMICRO)   :: ZWQ1, ZWQ2, ZWQ3, ZWQ4         !      for 
REAL,    DIMENSION(KMICRO,9) :: ZWQ5                           ! charge transfer
!
! variables used to select between common parameters between ICEx and LIMA
INTEGER :: IMOM_C, IMOM_R, IMOM_I, IMOM_S, IMOM_G, IMOM_H ! number of moments for each hydrometeor
INTEGER :: IGAMINC,              &
           IACCLBDAR, IACCLBDAS, &
           IDRYLBDAR, IDRYLBDAS, IDRYLBDAG
!
REAL    :: ZCEXVT,                                               &
           ZALPHAI, ZNUI, ZAI, ZBI, ZDS, ZDG, ZCXS, ZCXG,        &
           ZCOLIS, ZCOLEXIS, ZCOLIG, ZCOLEXIG, ZCOLSG, ZCOLEXSG, &
           ZACCINTP1S, ZACCINTP2S, ZACCINTP1R, ZACCINTP2R,       &
           ZDRYINTP1R, ZDRYINTP2R, ZDRYINTP1S, ZDRYINTP2S,       &
           ZDRYINTP1G, ZDRYINTP2G,                               &
           ZRIMINTP1,  ZRIMINTP2
REAL, DIMENSION(:), ALLOCATABLE :: ZRTMIN
!
! microphysical tendencies have to be transformed in 1D arrays
REAL, DIMENSION(KMICRO) :: ZRVHENI, ZRCHONI, ZRRHONG, ZRVDEPS, ZRIAGGS,      &
                           ZRIAUTS, ZRVDEPG, ZRCAUTR, ZRCACCR, ZRREVAV,      &
                           ZRIMLTC, ZRCBERI, ZRSMLTG, ZRRACCSS, ZRRACCSG,    &
                           ZRSACCRG, ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRICFRRG, &
                           ZRRCFRIG, ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG,     &
                           ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, ZRGMLTR
! optional microphysical tendencies
REAL, DIMENSION(:), ALLOCATABLE :: ZRCMLTSR, ZRICFRR, ZRVHENC, ZRCHINC, ZRVHONH,     &
                                   ZRRCVRC, ZRICNVI, ZRVDEPI, ZRSHMSI, ZRGHMGI,      &
                                   ZRICIBU, ZRIRDSF, ZRCCORR2, ZRRCORR2, ZRICORR2,   &
                                   ZRWETGH, ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH,      &
                                   ZRRWETH, ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH,      &
                                   ZRGDRYH, ZRHMLTR, ZRDRYHG 
!
!------------------------------------------------------------------
ASSOCIATE(XCEXVT_I=>ICED%XCEXVT, XRTMIN_I=>ICED%XRTMIN,                       &
                                XALPHAI_I=>ICED%XALPHAI, XNUI_I=>ICED%XNUI, XAI_I=>ICED%XAI, XBI_I=>ICED%XBI, &
                                XDS_I=>ICED%XDS, XDG_I=>ICED%XDG,                             &
                                XCXS_I=>ICED%XCXS, XCXG_I=>ICED%XCXG,                         &
                                XCOLIS_I=>ICEP%XCOLIS, XCOLEXIS_I=>ICEP%XCOLEXIS,             &
                                XCOLIG_I=>ICEP%XCOLIG, XCOLEXIG_I=>ICEP%XCOLEXIG,             &
                                XCOLSG_I=>ICEP%XCOLSG, XCOLEXSG_I=>ICEP%XCOLEXSG,             &
                                NGAMINC_I=>ICEP%NGAMINC,                                 &                                
                                NACCLBDAR_I=>ICEP%NACCLBDAR, NACCLBDAS_I=>ICEP%NACCLBDAS,     &
                                XACCINTP1S_I=>ICEP%XACCINTP1S, XACCINTP2S_I=>ICEP%XACCINTP2S, &
                                XACCINTP1R_I=>ICEP%XACCINTP1R, XACCINTP2R_I=>ICEP%XACCINTP2R, &
                                NDRYLBDAR_I=>ICEP%NDRYLBDAR, NDRYLBDAS_I=>ICEP%NDRYLBDAS,     &
                                NDRYLBDAG_I=>ICEP%NDRYLBDAG,                             &
                                XDRYINTP1R_I=>ICEP%XDRYINTP1R, XDRYINTP2R_I=>ICEP%XDRYINTP2R, &
                                XDRYINTP1S_I=>ICEP%XDRYINTP1S, XDRYINTP2S_I=>ICEP%XDRYINTP2S, &
                                XDRYINTP1G_I=>ICEP%XDRYINTP1G, XDRYINTP2G_I=>ICEP%XDRYINTP2G, &
                                XRIMINTP1_I=>ICEP%XRIMINTP1,   XRIMINTP2_I=>ICEP%XRIMINTP2    )
!
!*      1.    INITIALIZATIONS
!             ---------------
!
!*      1.1   compute the loop bounds
!
IIB = D%NIB
IIE = D%NIE
IJB = D%NJB
IJE = D%NJE
IKB = D%NKB
IKE = D%NKE
!
!
!*      1.2   select parameters between ICEx and LIMA
!
IF (HCLOUD(1:3) == 'ICE') THEN
  ZCEXVT = XCEXVT_I
  IMOM_C = 1
  IMOM_R = 1
  IMOM_I = 2 ! Ni is diagnostic and always available
  IMOM_S = 1
  IMOM_G = 1
  IF (KRR == 7) THEN
    IMOM_H = 1
  ELSE
    IMOM_H = 0
  END IF
ELSE IF (HCLOUD == 'LIMA') THEN
  ZCEXVT = XCEXVT_L
  IMOM_C = NMOM_C
  IMOM_R = NMOM_R
  IMOM_I = 2 ! Ni is diagnostic and always available
  IMOM_S = NMOM_S
  IMOM_G = NMOM_G
  IMOM_H = NMOM_H
END IF
!
ZRHO00 = CST%XP00 / (CST%XRD * PTHVREFZIKB)
ZCOR00 = ZRHO00**ZCEXVT
!
!!provisoire pour changement interface!!!IF (LINDUCTIVE) ALLOCATE (ZEFIELDW(KMICRO))
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.3   packing
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!  optimization by looking for locations where
!!provisoire pour changement interface!!!!  the microphysical fields are larger than a minimal value only !!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!IF (KMICRO >= 0) THEN
!!provisoire pour changement interface!!!  IMICRO = COUNTJV(ODMICRO(:,:,:), II1(:), II2(:), II3(:))
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! some microphysical tendencies are optional: the corresponding 1D arrays must be allocated
!!provisoire pour changement interface!!!  IF (HCLOUD(1:3) == 'ICE') THEN  ! ICE3 scheme
!!provisoire pour changement interface!!!    ALLOCATE(ZRCMLTSR(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRICFRR(IMICRO))
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN  ! LIMA scheme
!!provisoire pour changement interface!!!    ALLOCATE(ZRVHENC(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRCHINC(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRVHONH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRRCVRC(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRICNVI(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRVDEPI(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRSHMSI(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRGHMGI(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRICIBU(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRIRDSF(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRCCORR2(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRRCORR2(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRICORR2(IMICRO))
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!  IF (KRR == 7) THEN ! hail activated
!!provisoire pour changement interface!!!    ALLOCATE(ZRWETGH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRCWETH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRIWETH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRSWETH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRGWETH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRRWETH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRCDRYH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRRDRYH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRIDRYH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRSDRYH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRGDRYH(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRHMLTR(IMICRO))
!!provisoire pour changement interface!!!    ALLOCATE(ZRDRYHG(IMICRO))
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  DO JL = 1, IMICRO
!!provisoire pour changement interface!!!    ZZT(JL)      = PZT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRHODREF(JL) = PRHODREF(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRHOCOR(JL)  = (ZRHO00 / ZRHODREF(JL))**ZCEXVT
!!provisoire pour changement interface!!!    ZRHODJ(JL)   = PRHODJ(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZCIT(JL)     = PCIT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_C == 2) ZCCT(JL) = PCCT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_R == 2) ZCRT(JL) = PCRT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_S == 2) ZCST(JL) = PCST(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_G == 2) ZCGT(JL) = PCGT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_H == 2) ZCHT(JL) = PCHT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (IMOM_C == 1) ZCCT(JL) = 0.
!!provisoire pour changement interface!!!    IF (IMOM_R == 1) ZCRT(JL) = 0.
!!provisoire pour changement interface!!!    IF (IMOM_S == 1) ZCST(JL) = 0.
!!provisoire pour changement interface!!!    IF (IMOM_G == 1) ZCGT(JL) = 0.
!!provisoire pour changement interface!!!    IF (IMOM_H == 1) ZCHT(JL) = 0.
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZRVT(JL)  = PRVT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRCT(JL)  = PRCT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRRT(JL)  = PRRT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRIT(JL)  = PRIT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRST(JL)  = PRST(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZRGT(JL)  = PRGT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (KRR == 7) ZRHT(JL) = PRHT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQPIT(JL) = PQPIT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQNIT(JL) = PQNIT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQCT(JL)  = PQCT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQRT(JL)  = PQRT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQIT(JL)  = PQIT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQST(JL)  = PQST(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    ZQGT(JL)  = PQGT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    IF (KRR == 7) ZQHT(JL) = PQHT(II1(JL),II2(JL),II3(JL))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQPIS(JL) = PQPIS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQNIS(JL) = PQNIS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQCS(JL)  = PQCS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQRS(JL)  = PQRS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQIS(JL)  = PQIS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQSS(JL)  = PQSS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZQGS(JL)  = PQGS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    IF (KRR == 7) ZQHS(JL) = PQHS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IF (LINDUCTIVE) ZEFIELDW(JL) = XEFIELDW(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! microphysical tendencies
!!provisoire pour changement interface!!!    ZRVHENI(JL) = PRVHENI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRHONG(JL) = PRRHONG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRIMLTC(JL) = PRIMLTC(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCHONI(JL) = PRCHONI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRVDEPS(JL) = PRVDEPS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRIAGGS(JL) = PRIAGGS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRIAUTS(JL) = PRIAUTS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRVDEPG(JL) = PRVDEPG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCAUTR(JL) = PRCAUTR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCACCR(JL) = PRCACCR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRREVAV(JL) = PRREVAV(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCRIMSS(JL) = PRCRIMSS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCRIMSG(JL) = PRCRIMSG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRSRIMCG(JL) = PRSRIMCG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRACCSS(JL) = PRRACCSS(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRACCSG(JL) = PRRACCSG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRSACCRG(JL) = PRSACCRG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRSMLTG(JL) = PRSMLTG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRICFRRG(JL) = PRICFRRG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRCFRIG(JL) = PRRCFRIG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCWETG(JL) = PRCWETG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRIWETG(JL) = PRIWETG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRWETG(JL) = PRRWETG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRSWETG(JL) = PRSWETG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCDRYG(JL) = PRCDRYG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRIDRYG(JL) = PRIDRYG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRRDRYG(JL) = PRRDRYG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRSDRYG(JL) = PRSDRYG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRGMLTR(JL) = PRGMLTR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    ZRCBERI(JL) = PRCBERI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    IF (HCLOUD(1:3) == 'ICE') THEN
!!provisoire pour changement interface!!!      ZRCMLTSR(JL) = PRCMLTSR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRICFRR(JL)  = PRICFRR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!    IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!      ZCST(JL)    = PCST(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZCGT(JL)    = PCGT(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRVHENC(JL) = PRVHENC(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRCHINC(JL) = PRCHINC(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRVHONH(JL) = PRVHONH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRRCVRC(JL) = PRRCVRC(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRICNVI(JL) = PRICNVI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRVDEPI(JL) = PRVDEPI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRSHMSI(JL) = PRSHMSI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRGHMGI(JL) = PRGHMGI(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRICIBU(JL) = PRICIBU(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRIRDSF(JL) = PRIRDSF(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRCCORR2(JL) = PRCCORR2(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRRCORR2(JL) = PRRCORR2(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRICORR2(JL) = PRICORR2(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!    IF (KRR == 7) THEN
!!provisoire pour changement interface!!!      ZCHT(JL)    = PCHT(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRWETGH(JL) = PRWETGH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRCWETH(JL) = PRCWETH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRIWETH(JL) = PRIWETH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRSWETH(JL) = PRSWETH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRGWETH(JL) = PRGWETH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRRWETH(JL) = PRRWETH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRCDRYH(JL) = PRCDRYH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRRDRYH(JL) = PRRDRYH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRIDRYH(JL) = PRIDRYH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRSDRYH(JL) = PRSDRYH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRGDRYH(JL) = PRGDRYH(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRHMLTR(JL) = PRHMLTR(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!      ZRDRYHG(JL) = PRDRYHG(II1(JL), II2(JL), II3(JL))
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ZRHOCOR(:) = (ZRHO00 / ZRHODREF(:))**ZCEXVT
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.4   allocations for the non-inductive parameterizations
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'GARDI') THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZDELTALWC(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZFT(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'TAKAH' .OR.                              &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZEW(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'TEEWC') THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZDQ(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'TEEWC' )  THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSK(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIM(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIN(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSM(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSN(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'SAP98' .OR.                              &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZFQIAGGS(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZFQIDRYGBS(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZLBQSDRYGB1S(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZLBQSDRYGB2S(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZLBQSDRYGB3S(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'TERAR' .OR. &
!!provisoire pour changement interface!!!      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
!!provisoire pour changement interface!!!    ALLOCATE( ZRAR(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZDQ_IS(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZDQ_IG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZDQ_SG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIM_IS(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIN_IS(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIM_IG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNIN_IG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSK_SG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSM_SG(KMICRO) )
!!provisoire pour changement interface!!!    ALLOCATE( ZSAUNSN_SG(KMICRO) )
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.5   select parameters between ICEx and LIMA
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ALLOCATE(ZRTMIN(KRR))
!!provisoire pour changement interface!!!  IF (HCLOUD(1:3) == 'ICE') THEN
!!provisoire pour changement interface!!!! in ini_rain_ice, xrtmin is initialized with dimension 6 (hail not activated) or 7 (hail activated)
!!provisoire pour changement interface!!!    ZRTMIN(1:KRR) = XRTMIN_I(1:KRR)
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZALPHAI = XALPHAI_I
!!provisoire pour changement interface!!!    ZNUI = XNUI_I
!!provisoire pour changement interface!!!    ZAI = XAI_I
!!provisoire pour changement interface!!!    ZBI = XBI_I
!!provisoire pour changement interface!!!    ZDS = XDS_I
!!provisoire pour changement interface!!!    ZDG = XDG_I
!!provisoire pour changement interface!!!    ZCXS = XCXS_I
!!provisoire pour changement interface!!!    ZCXG = XCXG_I
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZCOLIS   = XCOLIS_I
!!provisoire pour changement interface!!!    ZCOLEXIS = XCOLEXIS_I
!!provisoire pour changement interface!!!    ZCOLIG   = XCOLIG_I
!!provisoire pour changement interface!!!    ZCOLEXIG = XCOLEXIG_I
!!provisoire pour changement interface!!!    ZCOLSG   = XCOLSG_I
!!provisoire pour changement interface!!!    ZCOLEXSG = XCOLEXSG_I
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IGAMINC = NGAMINC_I
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IACCLBDAR = NACCLBDAR_I
!!provisoire pour changement interface!!!    IACCLBDAS = NACCLBDAS_I
!!provisoire pour changement interface!!!    ZACCINTP1S = XACCINTP1S_I
!!provisoire pour changement interface!!!    ZACCINTP2S = XACCINTP2S_I
!!provisoire pour changement interface!!!    ZACCINTP1R = XACCINTP1R_I
!!provisoire pour changement interface!!!    ZACCINTP2R = XACCINTP2R_I
!!provisoire pour changement interface!!!    ! 
!!provisoire pour changement interface!!!    IDRYLBDAR = NDRYLBDAR_I
!!provisoire pour changement interface!!!    IDRYLBDAS = NDRYLBDAS_I
!!provisoire pour changement interface!!!    IDRYLBDAG = NDRYLBDAG_I
!!provisoire pour changement interface!!!    ZDRYINTP1R = XDRYINTP1R_I
!!provisoire pour changement interface!!!    ZDRYINTP2R = XDRYINTP2R_I
!!provisoire pour changement interface!!!    ZDRYINTP1S = XDRYINTP1S_I
!!provisoire pour changement interface!!!    ZDRYINTP2S = XDRYINTP2S_I
!!provisoire pour changement interface!!!    ZDRYINTP1G = XDRYINTP1G_I
!!provisoire pour changement interface!!!    ZDRYINTP2G = XDRYINTP2G_I
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZRIMINTP1 = XRIMINTP1_I
!!provisoire pour changement interface!!!    ZRIMINTP2 = XRIMINTP2_I
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!  ELSE IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!! in ini_lima, xrtmin is initialized with dimension 7
!!provisoire pour changement interface!!!    ZRTMIN(1:KRR) = XRTMIN_L(1:KRR)
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZALPHAI = XALPHAI_L
!!provisoire pour changement interface!!!    ZNUI = XNUI_L
!!provisoire pour changement interface!!!    ZAI = XAI_L
!!provisoire pour changement interface!!!    ZBI = XBI_L
!!provisoire pour changement interface!!!    ZDS = XDS_L
!!provisoire pour changement interface!!!    ZDG = XDG_L
!!provisoire pour changement interface!!!    ZCXS = XCXS_L
!!provisoire pour changement interface!!!    ZCXG = XCXG_L
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZCOLIS   = 0.25  ! variable not defined in LIMA, the value of ICEx is used 
!!provisoire pour changement interface!!!    ZCOLEXIS = XCOLEXIS_L
!!provisoire pour changement interface!!!    ZCOLIG   = XCOLIG_L
!!provisoire pour changement interface!!!    ZCOLEXIG = XCOLEXIG_L
!!provisoire pour changement interface!!!    ZCOLSG   = XCOLSG_L
!!provisoire pour changement interface!!!    ZCOLEXSG = XCOLEXSG_L
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IGAMINC = NGAMINC_L
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IACCLBDAR = NACCLBDAR_L
!!provisoire pour changement interface!!!    IACCLBDAS = NACCLBDAS_L
!!provisoire pour changement interface!!!    ZACCINTP1S = XACCINTP1S_L
!!provisoire pour changement interface!!!    ZACCINTP2S = XACCINTP2S_L
!!provisoire pour changement interface!!!    ZACCINTP1R = XACCINTP1R_L
!!provisoire pour changement interface!!!    ZACCINTP2R = XACCINTP2R_L
!!provisoire pour changement interface!!!    ! 
!!provisoire pour changement interface!!!    IDRYLBDAR = NDRYLBDAR_L
!!provisoire pour changement interface!!!    IDRYLBDAS = NDRYLBDAS_L
!!provisoire pour changement interface!!!    IDRYLBDAG = NDRYLBDAG_L
!!provisoire pour changement interface!!!    ZDRYINTP1R = XDRYINTP1R_L
!!provisoire pour changement interface!!!    ZDRYINTP2R = XDRYINTP2R_L
!!provisoire pour changement interface!!!    ZDRYINTP1S = XDRYINTP1S_L
!!provisoire pour changement interface!!!    ZDRYINTP2S = XDRYINTP2S_L
!!provisoire pour changement interface!!!    ZDRYINTP1G = XDRYINTP1G_L
!!provisoire pour changement interface!!!    ZDRYINTP2G = XDRYINTP2G_L
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZRIMINTP1 = XRIMINTP1_L
!!provisoire pour changement interface!!!    ZRIMINTP2 = XRIMINTP2_L
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.6   update the slope parameter of the distribution
!!provisoire pour changement interface!!!!*            and compute N_x if necessary
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD(1:3) == 'ICE') ZCCT(:) = 0.
!!provisoire pour changement interface!!!  CALL COMPUTE_LAMBDA(2, IMOM_C, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(2), ZRCT, ZCCT, ZLBDAC)
!!provisoire pour changement interface!!!  CALL COMPUTE_LAMBDA(3, IMOM_R, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(3), ZRRT, ZCRT, ZLBDAR)
!!provisoire pour changement interface!!!  CALL COMPUTE_LAMBDA(4, IMOM_I, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(4), ZRIT, ZCIT, ZLBDAI)
!!provisoire pour changement interface!!!  CALL COMPUTE_LAMBDA(5, IMOM_S, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(5), ZRST, ZCST, ZLBDAS)
!!provisoire pour changement interface!!!  CALL COMPUTE_LAMBDA(6, IMOM_G, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(6), ZRGT, ZCGT, ZLBDAG)
!!provisoire pour changement interface!!!  IF (KRR == 7) CALL COMPUTE_LAMBDA(7, IMOM_H, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(7), ZRHT, ZCHT, ZLBDAH)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.7   update the parameter e in the charge-diameter relationship
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! Compute e_x at time t
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(2, IMOM_C, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(2), ZRCT, ZQCT, ZECT, PLBDX=ZLBDAC, PCX=ZCCT)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(3, IMOM_R, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(3), ZRRT, ZQRT, ZERT, PLBDX=ZLBDAR, PCX=ZCRT)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(4, IMOM_I, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(4), ZRIT, ZQIT, ZEIT, PLBDX=ZLBDAI, PCX=ZCIT)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(5, IMOM_S, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(5), ZRST, ZQST, ZEST, PLBDX=ZLBDAS, PCX=ZCST)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(6, IMOM_G, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(6), ZRGT, ZQGT, ZEGT, PLBDX=ZLBDAG, PCX=ZCGT)
!!provisoire pour changement interface!!!    IF (KRR == 7) CALL ELEC_COMPUTE_EX(7, IMOM_H, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(7), ZRHT, ZQHT, ZEHT, PLBDX=ZLBDAH, PCX=ZCHT)
!!provisoire pour changement interface!!!  ELSE IF (HCLOUD(1:3) == 'ICE') THEN
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(2, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(2), ZRCT, ZQCT, ZECT)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(3, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(3), ZRRT, ZQRT, ZERT, PLBDX=ZLBDAR)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(4, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(4), ZRIT, ZQIT, ZEIT, PCX=ZCIT)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(5, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(5), ZRST, ZQST, ZEST, PLBDX=ZLBDAS)
!!provisoire pour changement interface!!!    CALL ELEC_COMPUTE_EX(6, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(6), ZRGT, ZQGT, ZEGT, PLBDX=ZLBDAG)
!!provisoire pour changement interface!!!    IF (KRR == 7) CALL ELEC_COMPUTE_EX(7, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(7), ZRHT, ZQHT, ZEHT, PLBDX=ZLBDAH)
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      1.8   initialization for the non-inductive charging process
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  SELECT CASE (CNI_CHARGING)
!!provisoire pour changement interface!!!    ! Initialization for the parameterization of Gardiner et al. (1995)
!!provisoire pour changement interface!!!    CASE ('GARDI')
!!provisoire pour changement interface!!!      CALL ELEC_INIT_NOIND_GARDI()
!!provisoire pour changement interface!!!      ! Save the effective water content
!!provisoire pour changement interface!!!      DO JL = 1, KMICRO
!!provisoire pour changement interface!!!        XEW(II1(JL),II2(JL),II3(JL)) = ZDELTALWC(JL) ! 
!!provisoire pour changement interface!!!      END DO      
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Initialization for the parameterizations of Saunders et al. (1991)
!!provisoire pour changement interface!!!    ! with and without anomalies, and Tsenova and Mitzeva (2009)
!!provisoire pour changement interface!!!    CASE ('SAUN1', 'SAUN2', 'TEEWC')
!!provisoire pour changement interface!!!      CALL ELEC_INIT_NOIND_EWC()
!!provisoire pour changement interface!!!      ! Save the effective water content
!!provisoire pour changement interface!!!      DO JL = 1, KMICRO
!!provisoire pour changement interface!!!        XEW(II1(JL),II2(JL),II3(JL)) = ZEW(JL) ! g/m3
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Initialization for the parameterizations of Saunders and Peck (1998), 
!!provisoire pour changement interface!!!    ! Brooks et al. (1997) and Tsenova and Mitzeva (2011)
!!provisoire pour changement interface!!!    CASE ('SAP98', 'BSMP1', 'BSMP2', 'TERAR')
!!provisoire pour changement interface!!!      CALL ELEC_INIT_NOIND_RAR()
!!provisoire pour changement interface!!!      ! Save the rime accretion rate (not recorded properly: 3 different RAR are computed !!!)
!!provisoire pour changement interface!!!      DO JL = 1, KMICRO
!!provisoire pour changement interface!!!        XEW(II1(JL),II2(JL),II3(JL)) = ZRAR(JL) ! g/m3
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Initialization for the parameterization of Takahashi (1978)
!!provisoire pour changement interface!!!    CASE ('TAKAH')
!!provisoire pour changement interface!!!      CALL ELEC_INIT_NOIND_TAKAH()
!!provisoire pour changement interface!!!      ! Save the effective water content
!!provisoire pour changement interface!!!      DO JL = 1, KMICRO
!!provisoire pour changement interface!!!        XEW(II1(JL),II2(JL),II3(JL)) = ZEW(JL) ! g/m3
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!  END SELECT
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.    COMPUTE THE SLOW COLD PROCESS SOURCES
!!provisoire pour changement interface!!!!             -------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.1   heterogeneous nucleation
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! --> rien n'est fait pour l'elec pour le moment
!!provisoire pour changement interface!!!! ICE3/4 : rvheni/rvhind
!!provisoire pour changement interface!!!! LIMA : rvhenc, rchinc, rvhonh
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.2   spontaneous freezing (rhong)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%INIT_PHY(D, 'SFR', &
!!provisoire pour changement interface!!!                            Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%INIT_PHY(D, 'SFR', &
!!provisoire pour changement interface!!!                            Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!  
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRRHONG(:) > 0. .AND. &
!!provisoire pour changement interface!!!         ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3))
!!provisoire pour changement interface!!!    ZWQ(:) = ZQRS(:)
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZQRS(:)
!!provisoire pour changement interface!!!    ZQRS(:) = 0.
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%END_PHY(D, 'SFR', &
!!provisoire pour changement interface!!!                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%END_PHY(D, 'SFR', &
!!provisoire pour changement interface!!!                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.3   cloud ice melting (rimltc)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%INIT_PHY(D, 'IMLT', &
!!provisoire pour changement interface!!!                            Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%INIT_PHY(D, 'IMLT', &
!!provisoire pour changement interface!!!                            Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  WHERE (ZRIMLTC(:) > 0.)
!!provisoire pour changement interface!!!    ZQCS(:) = ZQCS(:) + ZQIS(:)
!!provisoire pour changement interface!!!    ZQIS(:) = 0.
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%END_PHY(D, 'IMLT', &
!!provisoire pour changement interface!!!                           Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%END_PHY(D, 'IMLT', &
!!provisoire pour changement interface!!!                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.4   riming-conversion of the large sized aggregates into graupel ???
!!provisoire pour changement interface!!!! ancienne param => on calcule plutot cette tendance un peu plus loin ?
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.5   homogeneous nucleation (rchoni)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! CB : traitement different entre ice3 et lima --> a modifier eventuellement
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRCHONI(:) > 0. .AND.           &
!!provisoire pour changement interface!!!         ZRCT(:) > XRTMIN_ELEC(2) .AND.  &
!!provisoire pour changement interface!!!         ABS(ZQCT(:)) > XQTMIN(2) .AND. ABS(ZECT(:)) > ELECP%XECMIN)
!!provisoire pour changement interface!!!    ZWQ(:) = XQHON * ZECT(:) * ZRCHONI(:)
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQCS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQIS(:) = ZQIS(:) + ZWQ(:)
!!provisoire pour changement interface!!!    ZQCS(:) = ZQCS(:) - ZWQ(:)
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'HON', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'HON', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.6   deposition on snow/aggregates (rvdeps)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%INIT_PHY(D, 'DEPS', &
!!provisoire pour changement interface!!!                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%INIT_PHY(D, 'DEPS', &
!!provisoire pour changement interface!!!                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! Only the sublimation of snow/aggregates is considered (negative part of PRVDEPS)
!!provisoire pour changement interface!!!  WHERE (ZRVDEPS(:) < 0. .AND. & 
!!provisoire pour changement interface!!!         ZRST(:) > XRTMIN_ELEC(5) .AND. ABS(ZQST(:)) > XQTMIN(5))
!!provisoire pour changement interface!!!    ZWQ(:) = XCOEF_RQ_S * ZQST(:) * ZRVDEPS(:) / ZRST(:)
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ(:)) ),ZQSS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQSS(:)  = ZQSS(:)  - ZWQ(:)
!!provisoire pour changement interface!!!    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%END_PHY(D, 'DEPS', &
!!provisoire pour changement interface!!!                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%END_PHY(D, 'DEPS', &
!!provisoire pour changement interface!!!                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'DEPS', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.7   aggregation on snow/aggregates (riaggs)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRIAGGS, ZRIT, ZQIT, PTSTEP, &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_I, &
!!provisoire pour changement interface!!!                                ZWQ, ZQIS, ZQSS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'AGGS', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'AGGS', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.8   non-inductive charging during ice - snow collisions
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL ELEC_IAGGS_B()
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'NIIS', &
!!provisoire pour changement interface!!!                           Unpack( -zwq_ni(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'NIIS', &
!!provisoire pour changement interface!!!                           Unpack(  zwq_ni(:), mask = odmicro(:, :, :), field =  0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! Save the NI charging rate
!!provisoire pour changement interface!!!  DO JL = 1, KMICRO
!!provisoire pour changement interface!!!    XNI_IAGGS(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.9   autoconversion of r_i for r_s production (riauts/ricnvs)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRIAUTS, ZRIT, ZQIT, PTSTEP, &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_I, &
!!provisoire pour changement interface!!!                                ZWQ, ZQIS, ZQSS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'AUTS', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'AUTS', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.10  snow --> ice conversion (rscnvi)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRICNVI, ZRST, ZQST, PTSTEP, &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQSS, ZQIS)          
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'CNVI', &
!!provisoire pour changement interface!!!                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'CNVI', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.11  water vapor deposition on ice crystals (rvdepi)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%INIT_PHY(D, 'SUBI', &
!!provisoire pour changement interface!!!                              Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%INIT_PHY(D, 'SUBI', &
!!provisoire pour changement interface!!!                              Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZWQ(:) = 0.
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Only the sublimation of ice crystals is considered (negative part of PRVDEPI)
!!provisoire pour changement interface!!!    WHERE (ZRVDEPI(:) < 0. .AND. &
!!provisoire pour changement interface!!!           ZRIT(:) > XRTMIN_ELEC(4) .AND. ABS(ZQIT(:)) > XQTMIN(4))
!!provisoire pour changement interface!!!      ZWQ(:) = XCOEF_RQ_I * ZQIT(:) * ZRVDEPI(:) / ZRIT(:)
!!provisoire pour changement interface!!!      ZWQ(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQIS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQIS(:)  = ZQIS(:)  - ZWQ(:)
!!provisoire pour changement interface!!!      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%END_PHY(D, 'SUBI', &
!!provisoire pour changement interface!!!                             Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%END_PHY(D, 'SUBI', &
!!provisoire pour changement interface!!!                             Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )           
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'SUBI', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if    
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      2.12  water vapor deposition on graupel (rvdepg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%INIT_PHY(D, 'DEPG', &
!!provisoire pour changement interface!!!                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%INIT_PHY(D, 'DEPG', &
!!provisoire pour changement interface!!!                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! Only the sublimation of graupel is considered (negative part of PRVDEPG)
!!provisoire pour changement interface!!!  WHERE (ZRVDEPG(:) < 0. .AND. &
!!provisoire pour changement interface!!!         ZRGT(:) > XRTMIN_ELEC(6) .AND. ABS(ZQGT(:)) > XQTMIN(6))
!!provisoire pour changement interface!!!    ZWQ(:) = XCOEF_RQ_G * ZQGT(:) * ZRVDEPG(:) / ZRGT(:)
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQGS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:)  = ZQGS(:)  - ZWQ(:) 
!!provisoire pour changement interface!!!    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%END_PHY(D, 'DEPG', &
!!provisoire pour changement interface!!!                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%END_PHY(D, 'DEPG', &
!!provisoire pour changement interface!!!                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'DEPG', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      3.    COMPUTE THE WARM PROCESS SOURCES
!!provisoire pour changement interface!!!!             --------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      3.1   autoconversion of r_c for r_r production (rcautr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRCAUTR, ZRCT, ZQCT, PTSTEP,           &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
!!provisoire pour changement interface!!!                                ZWQ, ZQCS, ZQRS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'AUTO', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'AUTO', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      3.2   accretion of r_c for r_r production (rcaccr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRCACCR, ZRCT, ZQCT, PTSTEP,           &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
!!provisoire pour changement interface!!!                                ZWQ, ZQCS, ZQRS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'ACCR', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'ACCR', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      3.3   evaporation of raindrops (rrevav)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%INIT_PHY(D, 'REVA', &
!!provisoire pour changement interface!!!                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%INIT_PHY(D, 'REVA', &
!!provisoire pour changement interface!!!                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRREVAV(:) > 0. .AND. &
!!provisoire pour changement interface!!!         ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3))
!!provisoire pour changement interface!!!    ZWQ(:) = XCOEF_RQ_R * ZQRT(:) * ZRREVAV(:) / ZRRT(:)
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQRS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQRS(:)  = ZQRS(:)  - ZWQ(:)
!!provisoire pour changement interface!!!    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%END_PHY(D, 'REVA', &
!!provisoire pour changement interface!!!                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%END_PHY(D, 'REVA', &
!!provisoire pour changement interface!!!                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'REVA', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      3.4   conversion of drops to droplets (rrcvrc)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%INIT_PHY(D, 'R2C1', &
!!provisoire pour changement interface!!!                              Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%INIT_PHY(D, 'R2C1', &
!!provisoire pour changement interface!!!                              Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRRCVRC, ZRRT, ZQRT, PTSTEP,           &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(3), XQTMIN(3), XCOEF_RQ_R, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQRS, ZQCS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%END_PHY(D, 'R2C1', &
!!provisoire pour changement interface!!!                             Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%END_PHY(D, 'R2C1', &
!!provisoire pour changement interface!!!                             Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_s
!!provisoire pour changement interface!!!!             ---------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.1   cloud droplet riming of the aggregates
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%INIT_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%INIT_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqss(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%INIT_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.1.1 riming of the small sized aggregates (rcrimss)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRCRIMSS, ZRCT, ZQCT, PTSTEP,          &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
!!provisoire pour changement interface!!!                                ZWQ, ZQCS, ZQSS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.1.2 riming conversion of the large sized aggregates into graupel (rcrimsg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRCRIMSG(:) > 0. .AND.                                        &
!!provisoire pour changement interface!!!         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!         ABS(ZQCT(:)) > XQTMIN(2))
!!provisoire pour changement interface!!!    ZWQ(:) = XCOEF_RQ_C * ZQCT(:) * ZRCRIMSG(:) / ZRCT(:)      ! QCRIMSG
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQCS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ(:)
!!provisoire pour changement interface!!!    ZQCS(:) = ZQCS(:) - ZWQ(:)
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.1.3 riming conversion of the large sized aggregates into graupel (rsrimcg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  GMASK(:) = .FALSE.
!!provisoire pour changement interface!!!  IGMASK = 0
!!provisoire pour changement interface!!!  DO JL = 1, SIZE(GMASK)
!!provisoire pour changement interface!!!    IF (ZRSRIMCG(JL) > 0. .AND. ZZT(JL) < CST%XTT .AND. &
!!provisoire pour changement interface!!!        ZRCT(JL) > XRTMIN_ELEC(2) .AND. ZRST(JL) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!        ZLBDAS(JL) > 0.) THEN  !++cb-- 12/07/23 condition ajoutee pour eviter log(0)
!!provisoire pour changement interface!!!      IGMASK = IGMASK + 1
!!provisoire pour changement interface!!!      I1(IGMASK) = JL
!!provisoire pour changement interface!!!      GMASK(JL) = .TRUE.
!!provisoire pour changement interface!!!    ELSE
!!provisoire pour changement interface!!!      GMASK(JL) = .FALSE.
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ALLOCATE(ZVEC1(IGMASK))
!!provisoire pour changement interface!!!  ALLOCATE(ZVEC2(IGMASK))
!!provisoire pour changement interface!!!  ALLOCATE(IVEC2(IGMASK))
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! select the ZLBDAS
!!provisoire pour changement interface!!!  DO JL = 1, IGMASK
!!provisoire pour changement interface!!!    ZVEC1(JL) = ZLBDAS(I1(JL))
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  ! find the next lower indice for the ZLBDAS in the geometrical set of Lbda_s 
!!provisoire pour changement interface!!!  ! used to tabulate some moments of the incomplete gamma function
!!provisoire pour changement interface!!!  ZVEC2(1:IGMASK) = MAX( 1.00001, MIN( REAL(IGAMINC)-0.00001,           &
!!provisoire pour changement interface!!!                        ZRIMINTP1 * LOG( ZVEC1(1:IGMASK) ) + ZRIMINTP2 ) )
!!provisoire pour changement interface!!!  IVEC2(1:IGMASK) = INT( ZVEC2(1:IGMASK) )
!!provisoire pour changement interface!!!  ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL( IVEC2(1:IGMASK) )
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! perform the linear interpolation of the normalized "XFS"-moment of
!!provisoire pour changement interface!!!  ! the incomplete gamma function
!!provisoire pour changement interface!!!  ZVEC1(1:IGMASK) =  XGAMINC_RIM3( IVEC2(1:IGMASK)+1 ) *  ZVEC2(1:IGMASK)      &
!!provisoire pour changement interface!!!                   - XGAMINC_RIM3( IVEC2(1:IGMASK)   ) * (ZVEC2(1:IGMASK) - 1.0)
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  DO JL = 1, IGMASK
!!provisoire pour changement interface!!!    ZWQ(I1(JL)) = ZVEC1(JL)
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  DEALLOCATE(ZVEC1)
!!provisoire pour changement interface!!!  DEALLOCATE(ZVEC2)
!!provisoire pour changement interface!!!  DEALLOCATE(IVEC2)
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! riming-conversion of the large sized aggregates into graupeln (rsrimcg)
!!provisoire pour changement interface!!!  WHERE (ZRSRIMCG(:) > 0. .AND.                                        &
!!provisoire pour changement interface!!!         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!         ABS(ZQCT(:)) > XQTMIN(2) .AND. ABS(ZEST(:)) > XESMIN)
!!provisoire pour changement interface!!!    ZWQ(:) = XQSRIMCG * ZEST(:) * ZCST(:) *          &    ! QSRIMCG
!!provisoire pour changement interface!!!             ZLBDAS(:)**XEXQSRIMCG * (1. - ZWQ(:)) / &
!!provisoire pour changement interface!!!             (PTSTEP * ZRHODREF(:))
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ(:)) ),ZQSS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ(:)
!!provisoire pour changement interface!!!    ZQSS(:) = ZQSS(:) - ZWQ(:)
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%END_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqcs(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%END_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqss(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%END_PHY(D, 'RIM', &
!!provisoire pour changement interface!!!                           Unpack( zqgs(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.2   Hallett-Mossop ice multiplication process due to snow riming (rhmsi)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRSHMSI, ZRST, ZQST, PTSTEP,          &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_S, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQSS, ZQIS)       
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'HMS', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'HMS', &
!!provisoire pour changement interface!!!                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.3   Raindrop accretion onto the aggregates
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IGMASK = 0
!!provisoire pour changement interface!!!  DO JL = 1, SIZE(GMASK)
!!provisoire pour changement interface!!!    IF (ZRRT(JL) > ZRTMIN(3) .AND. ZLBDAR(JL) > 0. .AND. &
!!provisoire pour changement interface!!!        ZRST(JL) > ZRTMIN(5) .AND. ZLBDAS(JL) > 0.) THEN
!!provisoire pour changement interface!!!      IGMASK = IGMASK + 1
!!provisoire pour changement interface!!!      I1(IGMASK) = JL
!!provisoire pour changement interface!!!      GMASK(JL) = .TRUE.
!!provisoire pour changement interface!!!    ELSE
!!provisoire pour changement interface!!!      GMASK(JL) = .FALSE.
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  IF (IGMASK > 0) THEN
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ3(IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! select the (ZLBDAS,ZLBDAR) couplet
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZVEC1(JL) = ZLBDAS(I1(JL))
!!provisoire pour changement interface!!!      ZVEC2(JL) = ZLBDAR(I1(JL))
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! find the next lower indice for the ZLBDAS and for the ZLBDAR in the geometrical 
!!provisoire pour changement interface!!!    ! set of (Lbda_s,Lbda_r) couplet use to tabulate the kernels
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = MAX( 1.00001, MIN( REAL(IACCLBDAS)-0.00001,           &
!!provisoire pour changement interface!!!                          ZACCINTP1S * LOG( ZVEC1(1:IGMASK) ) + ZACCINTP2S ) )
!!provisoire pour changement interface!!!    IVEC1(1:IGMASK) = INT( ZVEC1(1:IGMASK) )
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL( IVEC1(1:IGMASK) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = MAX( 1.00001, MIN( REAL(IACCLBDAR)-0.00001,           &
!!provisoire pour changement interface!!!                          ZACCINTP1R * LOG( ZVEC2(1:IGMASK) ) + ZACCINTP2R ) )
!!provisoire pour changement interface!!!    IVEC2(1:IGMASK) = INT( ZVEC2(1:IGMASK) )
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL( IVEC2(1:IGMASK) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! perform the bilinear interpolation of the normalized kernels
!!provisoire pour changement interface!!!    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_RACCSS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
!!provisoire pour changement interface!!!    ZVECQ2(:) = BI_LIN_INTP_V(XKER_Q_RACCS,  IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
!!provisoire pour changement interface!!!    ZVECQ3(:) = BI_LIN_INTP_V(XKER_Q_SACCRG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
!!provisoire pour changement interface!!!    ZWQ1(:) = 0.
!!provisoire pour changement interface!!!    ZWQ2(:) = 0.
!!provisoire pour changement interface!!!    ZWQ3(:) = 0.
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZWQ1(I1(JL)) = ZVECQ1(JL)
!!provisoire pour changement interface!!!      ZWQ2(I1(JL)) = ZVECQ2(JL)
!!provisoire pour changement interface!!!      ZWQ3(I1(JL)) = ZVECQ3(JL)
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ1)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ2)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ3)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.3.1 raindrop accretion onto the small sized aggregates (rraccss)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ZWQ4(:)   = 0.
!!provisoire pour changement interface!!!    ZWQ5(:,:) = 0.
!!provisoire pour changement interface!!!    WHERE (ZRRACCSS(:) > 0. .AND. &
!!provisoire pour changement interface!!!           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!           ZCRT(:) > 0.             .AND. ZCST(:) > 0.             .AND. &
!!provisoire pour changement interface!!!           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.           .AND. &
!!provisoire pour changement interface!!!           ABS(ZERT(:)) > ELECP%XERMIN)  ! and zzt(:) < xtt ?
!!provisoire pour changement interface!!!      ZWQ4(:) = XFQRACCS * ZERT(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) * &
!!provisoire pour changement interface!!!                ZCRT(:) * ZCST(:)                                        * &
!!provisoire pour changement interface!!!               (XLBQRACCS1 * ZLBDAR(:)**(-2.0 - XFR)                     + &
!!provisoire pour changement interface!!!                XLBQRACCS2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAS(:)**(-1.0) + &
!!provisoire pour changement interface!!!                XLBQRACCS3 * ZLBDAR(:)**(-XFR)       * ZLBDAS(:)**(-2.0))
!!provisoire pour changement interface!!!      ZWQ5(:,1) = ZWQ4(:) * ZWQ1(:)                       ! QRACCSS
!!provisoire pour changement interface!!!      ZWQ5(:,1) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,1)
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) + ZWQ5(:,1)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.3.2 raindrop accretion-conversion of the large sized aggregates into graupel 
!!provisoire pour changement interface!!!!*            (rsaccrg & rraccsg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ZWQ5(:,2) = ZWQ2(:) * ZWQ4(:) ! QRACCS
!!provisoire pour changement interface!!!    WHERE (ZRRACCSG(:) > 0. .AND. &
!!provisoire pour changement interface!!!           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.)
!!provisoire pour changement interface!!!      ZWQ5(:,3) = ZWQ5(:,2) - ZWQ5(:,1)         ! QRACCSG
!!provisoire pour changement interface!!!      ZWQ5(:,3) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,3)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,3)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    WHERE (ZRSACCRG(:) > 0. .AND.                                        &
!!provisoire pour changement interface!!!           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
!!provisoire pour changement interface!!!           ZCRT(:) > 0.             .AND. ZCST(:) > 0.             .AND. &
!!provisoire pour changement interface!!!           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.           .AND. &
!!provisoire pour changement interface!!!           ABS(ZEST(:)) > XESMIN)
!!provisoire pour changement interface!!!      ZWQ5(:,4) = ZWQ3(:) * XFQRACCS * ZEST(:) *                               &
!!provisoire pour changement interface!!!                  ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *                        &
!!provisoire pour changement interface!!!                  ZCRT(:) * ZCST(:) *                                          &
!!provisoire pour changement interface!!!                 (XLBQSACCRG1 * ZLBDAS(:)**(-2.0 - XFS) +                      &
!!provisoire pour changement interface!!!                  XLBQSACCRG2 * ZLBDAS(:)**(-1.0 - XFS) * ZLBDAR(:)**(-1.0) +  &
!!provisoire pour changement interface!!!                  XLBQSACCRG3 * ZLBDAS(:)**(-XFS)       * ZLBDAR(:)**(-2.0)) ! QSACCR
!!provisoire pour changement interface!!!      ZWQ5(:,4) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQSS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,4)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,4)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'ACC', &
!!provisoire pour changement interface!!!                             Unpack( (-zwq5(:,1) - zwq5(:,3)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'ACC', &
!!provisoire pour changement interface!!!                             Unpack( ( zwq5(:,1) - zwq5(:,4)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'ACC', &
!!provisoire pour changement interface!!!                             Unpack( ( zwq5(:,3) + zwq5(:,4)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if    
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!  END IF  ! end if igmask>0
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.4   conversion-melting of the aggregates (rsmltg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRSMLTG, ZRST, ZQST, PTSTEP,           &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
!!provisoire pour changement interface!!!                                ZWQ, ZQSS, ZQGS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'CMEL', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'CMEL', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      4.5   cloud droplet collection onto aggregates by positive temperature (rcmltsr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD(1:3) == 'ICE') THEN
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRCMLTSR, ZRCT, ZQCT, PTSTEP,           &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQCS, ZQRS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'CMEL', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'CMEL', &
!!provisoire pour changement interface!!!                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_g
!!provisoire pour changement interface!!!!             ---------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.1   rain contact freezing (ricfrrg, rrcfrig, ricfrr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%INIT_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%INIT_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%INIT_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRRCFRIG(:) > 0. .AND.                                        &
!!provisoire pour changement interface!!!         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
!!provisoire pour changement interface!!!         ZCRT(:) > 0.             .AND.                                &
!!provisoire pour changement interface!!!         ABS(ZERT(:)) > ELECP%XERMIN    .AND. ABS(ZQRT(:)) > XQTMIN(3))
!!provisoire pour changement interface!!!    ZWQ(:) = XQRCFRIG * ZLBDAR(:)**XEXQRCFRIG * ZCIT(:) * ZCRT(:) * &
!!provisoire pour changement interface!!!             ZERT(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))  ! QRCFRIG
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQRS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ(:)
!!provisoire pour changement interface!!!    ZQRS(:) = ZQRS(:) - ZWQ(:)
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ZWQ(:) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRICFRRG(:) > 0. .AND.                                        &
!!provisoire pour changement interface!!!         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
!!provisoire pour changement interface!!!         ABS(ZQIT(:)) > XQTMIN(4))
!!provisoire pour changement interface!!!    ZWQ(:) = XCOEF_RQ_I * ZQIT(:) * ZRICFRRG(:) / ZRIT(:)     ! QICFRRG
!!provisoire pour changement interface!!!    ZWQ(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQIS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ(:)
!!provisoire pour changement interface!!!    ZQIS(:) = ZQIS(:) - ZWQ(:)
!!provisoire pour changement interface!!!  ENDWHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!++CB-- 16/06/2022 il manque le traitement de qricfrr
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%END_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%END_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%END_PHY(D, 'CFRZ', &
!!provisoire pour changement interface!!!                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2   graupel dry growth (qcdryg, qrdryg, qidryg & qsdryg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ5(:,:) = 0.
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.1 compute qcdryg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  WHERE (ZRCDRYG(:) > 0. .AND.                                         &
!!provisoire pour changement interface!!!         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
!!provisoire pour changement interface!!!         ABS(ZQCT(:)) > XQTMIN(2)) 
!!provisoire pour changement interface!!!    ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCDRYG(:) / ZRCT(:)
!!provisoire pour changement interface!!!    ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ5(:,1)
!!provisoire pour changement interface!!!  ENDWHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.2 compute qidryg = qidryg_coal + qidryg_boun
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  WHERE (ZRIDRYG(:) > 0. .AND.                                         &
!!provisoire pour changement interface!!!         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
!!provisoire pour changement interface!!!         ABS(ZQIT(:)) > XQTMIN(4))
!!provisoire pour changement interface!!!    ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIDRYG(:) / ZRIT(:)         ! QIDRYG_coal
!!provisoire pour changement interface!!!    ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
!!provisoire pour changement interface!!!    ZQGS(:) = ZQGS(:) + ZWQ5(:,2)
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.3 compute non-inductive charging durig ice - graupel collisions
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ! charge separation during collision between ice and graupel
!!provisoire pour changement interface!!!  CALL ELEC_IDRYG_B()  ! QIDRYG_boun
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'NIIG', &
!!provisoire pour changement interface!!!                            Unpack( -zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'NIIG', &
!!provisoire pour changement interface!!!                            Unpack(  zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! Save the NI charging rate
!!provisoire pour changement interface!!!  DO JL = 1, KMICRO
!!provisoire pour changement interface!!!    XNI_IDRYG(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.4 compute qsdryg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IGMASK = 0
!!provisoire pour changement interface!!!  DO JL = 1, SIZE(GMASK)
!!provisoire pour changement interface!!!    IF (ZRST(JL) > ZRTMIN(5) .AND. ZLBDAS(JL) > 0. .AND. &
!!provisoire pour changement interface!!!        ZRGT(JL) > ZRTMIN(6) .AND. ZLBDAG(JL) > 0.) THEN
!!provisoire pour changement interface!!!      IGMASK = IGMASK + 1
!!provisoire pour changement interface!!!      I1(IGMASK) = JL
!!provisoire pour changement interface!!!      GMASK(JL) = .TRUE.
!!provisoire pour changement interface!!!    ELSE
!!provisoire pour changement interface!!!      GMASK(JL) = .FALSE.
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  IF (IGMASK > 0) THEN
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ3(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ4(IGMASK))
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! select the (ZLBDAG,ZLBDAS) couplet
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZVEC1(JL) = ZLBDAG(I1(JL))
!!provisoire pour changement interface!!!      ZVEC2(JL) = ZLBDAS(I1(JL))
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  ! find the next lower indice for the ZLBDAG and for the ZLBDAS in the geometrical set 
!!provisoire pour changement interface!!!  ! of (Lbda_g,Lbda_s) couplet use to tabulate the SDRYG-kernel
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = MAX(1.00001, MIN(REAL(IDRYLBDAG)-0.00001,  &
!!provisoire pour changement interface!!!                          ZDRYINTP1G*LOG(ZVEC1(1:IGMASK))+ZDRYINTP2G))
!!provisoire pour changement interface!!!    IVEC1(1:IGMASK) = INT(ZVEC1(1:IGMASK) )
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL(IVEC1(1:IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAS)-0.00001, &
!!provisoire pour changement interface!!!                          ZDRYINTP1S*LOG(ZVEC2(1:IGMASK))+ZDRYINTP2S))
!!provisoire pour changement interface!!!    IVEC2(1:IGMASK) = INT(ZVEC2(1:IGMASK))
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL(IVEC2(1:IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! perform the bilinear interpolation of the normalized QSDRYG-kernels
!!provisoire pour changement interface!!!    ! normalized Q-SDRYG-kernel 
!!provisoire pour changement interface!!!    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_SDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
!!provisoire pour changement interface!!!    ZWQ5(:,3) = 0.  ! normalement pas utile
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZWQ5(I1(JL),3) = ZVECQ1(JL)
!!provisoire pour changement interface!!!    END DO    
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! normalized Q-???-kernel
!!provisoire pour changement interface!!!    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAUN1' .OR. &
!!provisoire pour changement interface!!!        CNI_CHARGING == 'SAUN2' .OR. CNI_CHARGING == 'SAP98' .OR. &
!!provisoire pour changement interface!!!        CNI_CHARGING == 'GARDI' .OR.                              &
!!provisoire pour changement interface!!!        CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
!!provisoire pour changement interface!!!        CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
!!provisoire pour changement interface!!!      ZVECQ2(:)  = BI_LIN_INTP_V(XKER_Q_LIMSG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK) 
!!provisoire pour changement interface!!!      ZWQ5(:,4) = 0.  ! normalement pas utile
!!provisoire pour changement interface!!!      DO JL = 1, IGMASK
!!provisoire pour changement interface!!!        ZWQ5(I1(JL),4) = ZVECQ2(JL)
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! normalized Q-SDRYG-bouncing kernel
!!provisoire pour changement interface!!!    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'HELFA' .OR. &
!!provisoire pour changement interface!!!        CNI_CHARGING == 'GARDI') THEN
!!provisoire pour changement interface!!!      ZVECQ3(:)  = BI_LIN_INTP_V(XKER_Q_SDRYGB,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
!!provisoire pour changement interface!!!      ZWQ5(:,5) = 0.  ! normalement pas utile
!!provisoire pour changement interface!!!      DO JL = 1, IGMASK
!!provisoire pour changement interface!!!        ZWQ5(I1(JL),5) = ZVECQ3(JL)
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    ELSE
!!provisoire pour changement interface!!!      ZVECQ3(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB1,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
!!provisoire pour changement interface!!!      ZVECQ4(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB2,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
!!provisoire pour changement interface!!!      ZWQ5(:,6:7) = 0.  ! normalement pas utile
!!provisoire pour changement interface!!!      DO JL = 1, IGMASK
!!provisoire pour changement interface!!!        ZWQ5(I1(JL),6) = ZVECQ3(JL) ! Dvqsgmn if charge>0
!!provisoire pour changement interface!!!        ZWQ5(I1(JL),7) = ZVECQ4(JL) ! Dvqsgmn if charge<0
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    ENDIF
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ1)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ2)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ3)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ4)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!++CB-- CALCULER E_SG ICI POUR EVITER DES CALCULS REDONDANTS
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! compute QSDRYG_coal
!!provisoire pour changement interface!!!    WHERE (ZRSDRYG(:) > 0 .AND.                                & !GDRY(:) .AND. &
!!provisoire pour changement interface!!!           ZRST(:) > XRTMIN_ELEC(5) .AND.                      &
!!provisoire pour changement interface!!!           ZLBDAS(:) > 0.           .AND. ZLBDAG(:) > 0. .AND. &
!!provisoire pour changement interface!!!           ABS(ZQST(:)) > XQTMIN(5) .AND. ABS(ZEST(:)) > XESMIN)
!!provisoire pour changement interface!!!      ZWQ5(:,3) = ZWQ5(:,3) * XFQSDRYG *                                   &
!!provisoire pour changement interface!!!                  ZCOLSG * EXP(ZCOLEXSG * (ZZT(:) - CST%XTT)) *                &
!!provisoire pour changement interface!!!                  ZEST(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *          & 
!!provisoire pour changement interface!!!                  ZCGT(:) * ZCST(:) *                                      &
!!provisoire pour changement interface!!!                 (XLBQSDRYG1 * ZLBDAS(:)**(-2.0-XFS) +                     &
!!provisoire pour changement interface!!!                  XLBQSDRYG2 * ZLBDAS(:)**(-1.0-XFS) * ZLBDAG(:)**(-1.0) + &
!!provisoire pour changement interface!!!                  XLBQSDRYG3 * ZLBDAS(:)**(-XFS)     * ZLBDAG(:)**(-2.0)) ! QSDRYG_coal
!!provisoire pour changement interface!!!      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,3)
!!provisoire pour changement interface!!!    ELSEWHERE
!!provisoire pour changement interface!!!      ZWQ5(:,3) = 0.
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.5 compute non-inductive charging during snow - graupel collisions
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ! compute QSDRYG_boun
!!provisoire pour changement interface!!!    CALL ELEC_SDRYG_B()
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'NISG', &
!!provisoire pour changement interface!!!                              Unpack( -zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'NISG', &
!!provisoire pour changement interface!!!                              Unpack(  zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Save the NI charging rate
!!provisoire pour changement interface!!!    DO JL = 1, KMICRO
!!provisoire pour changement interface!!!      XNI_SDRYG(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!  END IF  ! end if igmask>0
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.2.6 compute qrdryg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IGMASK = 0
!!provisoire pour changement interface!!!  GMASK(:) = .FALSE.
!!provisoire pour changement interface!!!  DO JL = 1, SIZE(GMASK)
!!provisoire pour changement interface!!!    IF (ZRRT(JL) > ZRTMIN(3) .AND. ZLBDAR(JL) > 0. .AND. &
!!provisoire pour changement interface!!!        ZRGT(JL) > ZRTMIN(6) .AND. ZLBDAG(JL) > 0.) THEN
!!provisoire pour changement interface!!!      IGMASK = IGMASK + 1
!!provisoire pour changement interface!!!      I1(IGMASK) = JL
!!provisoire pour changement interface!!!      GMASK(JL) = .TRUE.
!!provisoire pour changement interface!!!    ELSE
!!provisoire pour changement interface!!!      GMASK(JL) = .FALSE.
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  !
!!provisoire pour changement interface!!!  IF (IGMASK > 0) THEN
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC1(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(IVEC2(IGMASK))
!!provisoire pour changement interface!!!    ALLOCATE(ZVECQ1(IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! select the (ZLBDAG,ZLBDAR) couplet
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZVEC1(JL) = ZLBDAG(I1(JL))
!!provisoire pour changement interface!!!      ZVEC2(JL) = ZLBDAR(I1(JL))
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! find the next lower indice for the ZLBDAG and for the ZLBDAR in the geometrical set 
!!provisoire pour changement interface!!!    ! of (Lbda_g,Lbda_r) couplet use to tabulate the QDRYG-kernel
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAG)-0.00001,           &
!!provisoire pour changement interface!!!                         ZDRYINTP1G*LOG(ZVEC1(1:IGMASK))+ZDRYINTP2G))
!!provisoire pour changement interface!!!    IVEC1(1:IGMASK) = INT(ZVEC1(1:IGMASK))
!!provisoire pour changement interface!!!    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL(IVEC1(1:IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAR)-0.00001,           &
!!provisoire pour changement interface!!!                          ZDRYINTP1R*LOG(ZVEC2(1:IGMASK))+ZDRYINTP2R))
!!provisoire pour changement interface!!!    IVEC2(1:IGMASK) = INT(ZVEC2(1:IGMASK))
!!provisoire pour changement interface!!!    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL(IVEC2(1:IGMASK))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! perform the bilinear interpolation of the normalized RDRYG-kernel
!!provisoire pour changement interface!!!    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_RDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
!!provisoire pour changement interface!!!    ZWQ5(:,4) = 0.
!!provisoire pour changement interface!!!    DO JL = 1, IGMASK
!!provisoire pour changement interface!!!      ZWQ5(I1(JL),4) = ZVECQ1(JL)
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC1)
!!provisoire pour changement interface!!!    DEALLOCATE(IVEC2)
!!provisoire pour changement interface!!!    DEALLOCATE(ZVECQ1)  
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! compute QRDRYG
!!provisoire pour changement interface!!!    WHERE (ZRRDRYG(:) > 0. .AND.                                         &
!!provisoire pour changement interface!!!           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. & 
!!provisoire pour changement interface!!!           ZCRT(:) > 0.             .AND. ZCGT(:) > 0.             .AND. & 
!!provisoire pour changement interface!!!           ZLBDAR(:) > 0.           .AND. ZLBDAG(:) > 0.           .AND. &
!!provisoire pour changement interface!!!           ABS(ZERT(:)) > ELECP%XERMIN .AND. ABS(ZQRT(:)) > XQTMIN(3))
!!provisoire pour changement interface!!!      ZWQ5(:,4) = ZWQ5(:,4) * XFQRDRYG *                                       &
!!provisoire pour changement interface!!!                  ZRHODREF(:)**(-ZCEXVT) *                                     &
!!provisoire pour changement interface!!!                  ZERT(:) * ZCGT(:) * ZCRT(:) *                                &
!!provisoire pour changement interface!!!                 (XLBQRDRYG1 * ZLBDAR(:)**(-2.0 - XFR) +                       &
!!provisoire pour changement interface!!!                  XLBQRDRYG2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAG(:)**(-1.0) +   &
!!provisoire pour changement interface!!!                  XLBQRDRYG3 * ZLBDAR(:)**(-XFR)       * ZLBDAG(:)**(-2.0))
!!provisoire pour changement interface!!!      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,4)
!!provisoire pour changement interface!!!    ELSEWHERE
!!provisoire pour changement interface!!!      ZWQ5(:,4) = 0.
!!provisoire pour changement interface!!!    ENDWHERE
!!provisoire pour changement interface!!!!    ZRDRYG(:) = ZWQ5(:,1) + ZWQ5(:,2) + ZWQ5(:,3) + ZWQ5(:,4)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'DRYG', &
!!provisoire pour changement interface!!!                              Unpack( -zwq5(:,1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'DRYG', &
!!provisoire pour changement interface!!!                              Unpack( -zwq5(:,4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'DRYG', &
!!provisoire pour changement interface!!!                              Unpack( -zwq5(:,2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'DRYG', &
!!provisoire pour changement interface!!!                              Unpack( -zwq5(:,3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'DRYG', &
!!provisoire pour changement interface!!!                              Unpack( (zwq5(:,1) + zwq5(:,2) + zwq5(:,3) + zwq5(:,4)) * zrhodj(:), &
!!provisoire pour changement interface!!!                              mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  END IF  ! end if igmask>0
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.3   Hallett-Mossop ice multiplication process due to graupel riming (rhmgi)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRGHMGI, ZRGT, ZQGT, PTSTEP,          &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(6), XQTMIN(6), XCOEF_RQ_G, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQGS, ZQIS)          
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'HMG', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'HMG', &
!!provisoire pour changement interface!!!                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4   graupel wet growth (rcwetg, rrwetg, riwetg & rswetg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4.1 compute qcwetg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ5(:,5) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRCWETG(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2) .AND. ABS(ZQCT(:)) > XQTMIN(2) .AND. &
!!provisoire pour changement interface!!!         ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!    ZWQ5(:,5) = XCOEF_RQ_C * ZRCWETG(:) * ZQCT(:) / ZRCT(:)
!!provisoire pour changement interface!!!    ZWQ5(:,5) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQCS(:) )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4.1 compute qiwetg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ5(:,6) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRIWETG(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4) .AND. ABS(ZQIT(:)) > XQTMIN(4) .AND. &
!!provisoire pour changement interface!!!         ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!    ZWQ5(:,6) = XCOEF_RQ_I * ZRIWETG(:) * ZQIT(:) / ZRIT(:)
!!provisoire pour changement interface!!!    ZWQ5(:,6) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,6)) ),ZQIS(:) )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4.2 compute qswetg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ5(:,7) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRSWETG(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. ABS(ZQST(:)) > XQTMIN(5) .AND. &
!!provisoire pour changement interface!!!         ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!    ZWQ5(:,7) = XCOEF_RQ_S * ZRSWETG(:) * ZQST(:) / ZRST(:)
!!provisoire pour changement interface!!!    ZWQ5(:,7) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,7)) ),ZQSS(:) )
!!provisoire pour changement interface!!!  END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4.3 compute qrwetg
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  ZWQ5(:,8) = 0.
!!provisoire pour changement interface!!!  WHERE (ZRRWETG(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3) .AND. &
!!provisoire pour changement interface!!!         ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!    ZWQ5(:,8) = XCOEF_RQ_R * ZQRT(:) * ZRRWETG(:) / ZRRT(:)
!!provisoire pour changement interface!!!    ZWQ5(:,8) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,8)) ),ZQRS(:) )
!!provisoire pour changement interface!!!  ENDWHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.4.4 conversion of graupel into hail (rwetgh)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%INIT_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                            Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (KRR == 7) THEN
!!provisoire pour changement interface!!!    ZWQ5(:,9) = 0.
!!provisoire pour changement interface!!!    WHERE (ZRWETGH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. ABS(ZQGT(:)) > XQTMIN(6))
!!provisoire pour changement interface!!!      ZWQ5(:,9) = XCOEF_RQ_G * ZQGT(:) * ZRWETGH(:) / ZRGT(:)
!!provisoire pour changement interface!!!      ZWQ5(:,9) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,9)) ),ZQGS(:) )
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRCWETG(:) > 0. .OR. ZRRWETG(:) > 0. .OR. ZRIWETG(:) > 0. .OR. &
!!provisoire pour changement interface!!!           ZRSWETG(:) > 0. .OR. ZRWETGH(:) > 0.)
!!provisoire pour changement interface!!!      ZQCS(:) = ZQCS(:) - ZWQ5(:,5)
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,8)
!!provisoire pour changement interface!!!      ZQIS(:) = ZQIS(:) - ZWQ5(:,6)
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,7)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,5) + ZWQ5(:,8) + ZWQ5(:,6) + ZWQ5(:,7) - ZWQ5(:,9)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,9)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!  ELSE IF (KRR == 6) THEN
!!provisoire pour changement interface!!!    WHERE (ZRCWETG(:) > 0. .OR. ZRRWETG(:) > 0. .OR. ZRIWETG(:) > 0. .OR. &
!!provisoire pour changement interface!!!           ZRSWETG(:) > 0.)
!!provisoire pour changement interface!!!      ZQCS(:) = ZQCS(:) - ZWQ5(:,5)
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,8)
!!provisoire pour changement interface!!!      ZQIS(:) = ZQIS(:) - ZWQ5(:,6)
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,7)
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) + ZWQ5(:,5) + ZWQ5(:,8) + ZWQ5(:,6) + ZWQ5(:,7)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack( -zwq5(:,5) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack( -zwq5(:,8) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack( -zwq5(:,6) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack( -zwq5(:,7) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%END_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack(  zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    if ( krr == 7 ) &
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 )%PTR%ADD_PHY(D, 'WETG', &
!!provisoire pour changement interface!!!                           Unpack( zwq5(:,9) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.5   compute charge separation by the inductive mechanism
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! Computation of the charge transfer rate during inductive mechanism
!!provisoire pour changement interface!!!! Only the bouncing droplet-graupel collision when the graupel is in the dry 
!!provisoire pour changement interface!!!! growth mode is considered
!!provisoire pour changement interface!!!! The electric field is limited to 100 kV/m
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (LINDUCTIVE) THEN
!!provisoire pour changement interface!!!    ZWQ(:) = 0.
!!provisoire pour changement interface!!!    GMASK(:) = ZRCDRYG(:) > 0. 
!!provisoire pour changement interface!!!    IGMASK   = COUNT(GMASK(:))
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    IF (IGMASK > 0) THEN
!!provisoire pour changement interface!!!      ZWQ(:) = 0.
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      WHERE (GMASK(:) .AND.                                      &
!!provisoire pour changement interface!!!             ZEFIELDW(:) /= 0. .AND. ABS(ZEGT(:)) > ELECP%XEGMIN .AND. &
!!provisoire pour changement interface!!!             ZLBDAG(:) > 0. .AND. ZCGT(:) > 0. .AND.             &
!!provisoire pour changement interface!!!             ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRCT(:) > XRTMIN_ELEC(2))
!!provisoire pour changement interface!!!        ZWQ(:) = XIND1 * ZCGT(:) * ZRHOCOR(:) *                             &
!!provisoire pour changement interface!!!                (XIND2 * SIGN(MIN(100.E3, ABS(ZEFIELDW(:))), ZEFIELDW(:)) * &
!!provisoire pour changement interface!!!                 ZLBDAG(:) **(-2.-ZDG) -                                    &
!!provisoire pour changement interface!!!                 XIND3 * ZEGT(:) * ZLBDAG(:)**(-XFG-ZDG))
!!provisoire pour changement interface!!!        ZWQ(:) = ZWQ(:) / ZRHODREF(:)
!!provisoire pour changement interface!!!        !
!!provisoire pour changement interface!!!        ZQGS(:) = ZQGS(:) + ZWQ(:)
!!provisoire pour changement interface!!!        ZQCS(:) = ZQCS(:) - ZWQ(:)
!!provisoire pour changement interface!!!      END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!      if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!        CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'INCG', &
!!provisoire pour changement interface!!!                               Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!        CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'INCG', &
!!provisoire pour changement interface!!!                               Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      end if    
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ! Save the inductive charging rate
!!provisoire pour changement interface!!!      DO JL = 1, KMICRO
!!provisoire pour changement interface!!!        XIND_RATE(II1(JL),II2(JL),II3(JL)) = ZWQ(JL) * ZRHODREF(JL) ! C/m3/s
!!provisoire pour changement interface!!!      END DO
!!provisoire pour changement interface!!!    END IF
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ! Save the inductive charging rate
!!provisoire pour changement interface!!!    DO JL = 1, KMICRO
!!provisoire pour changement interface!!!      XIND_RATE(II1(JL),II2(JL),II3(JL)) = ZWQ(JL) * ZRHODREF(JL) ! C/m3/s
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      5.6   melting of the graupel (rgmltr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRGMLTR, ZRGT, ZQGT, PTSTEP, &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(6), XQTMIN(6), XCOEF_RQ_G, &
!!provisoire pour changement interface!!!                                ZWQ, ZQGS, ZQRS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'GMLT', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'GMLT', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      6.    COMPUTE THE OPTIONAL SECONDARY ICE PRODUCTION
!!provisoire pour changement interface!!!!             ---------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! dans un premier temps, on considere que la charge echangee est proportionnelle
!!provisoire pour changement interface!!!! a la masse echangee
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      6.1   collisional ice breakup (cibu)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA' .AND. LCIBU) &
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRICIBU, ZRST, ZQST, PTSTEP,           &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQSS, ZQIS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      6.2   raindrop shattering freezing (rdsf)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA' .AND. LRDSF) &
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRIRDSF, ZRRT, ZQRT, PTSTEP,           &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(3), XQTMIN(3), XCOEF_RQ_R, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQRS, ZQIS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      7.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_h
!!provisoire pour changement interface!!!!             ---------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (KRR == 7) THEN
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      7.1   wet growth of hail (qcweth, qrweth, qiweth, qsweth, qgweth)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ZWQ5(:,:) = 0.
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRCWETH(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2))
!!provisoire pour changement interface!!!      ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCWETH(:) / ZRCT(:)        ! QCWETH
!!provisoire pour changement interface!!!      ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,1)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRIWETH(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4))
!!provisoire pour changement interface!!!      ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIWETH(:) / ZRIT(:)        ! QIWETH
!!provisoire pour changement interface!!!      ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,2)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRSWETH(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5))
!!provisoire pour changement interface!!!      ZWQ5(:,3) = XCOEF_RQ_S * ZQST(:) * ZRSWETH(:) / ZRST(:)        ! QSWETH
!!provisoire pour changement interface!!!      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,3)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRGWETH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!      ZWQ5(:,5) = XCOEF_RQ_G * ZQGT(:) * ZRGWETH(:) / ZRGT(:)        ! QGWETH
!!provisoire pour changement interface!!!      ZWQ5(:,5) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQGS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) - ZWQ5(:,5)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,5)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRRWETH(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3))
!!provisoire pour changement interface!!!      ZWQ5(:,4) = XCOEF_RQ_R * ZQRT(:) * ZRRWETH(:) / ZRRT(:)        ! QRWETH
!!provisoire pour changement interface!!!      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,4)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'WETH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'WETH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'WETH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'WETH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'WETH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 5) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 )%PTR%ADD_PHY(D, 'WETH',                      &
!!provisoire pour changement interface!!!                             Unpack( ( zwq5(:, 1) + zwq5(:, 2) + zwq5(:, 3) + zwq5(:, 4) + zwq5(:, 5) ) &
!!provisoire pour changement interface!!!                                       * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      7.2   dry growth of hail (qcdryh, qrdryh, qidryh, qsdryh, qgdryh)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ZWQ5(:,:) = 0.
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRCDRYH(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2))
!!provisoire pour changement interface!!!      ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCDRYH(:) / ZRCT(:)        ! QCDRYH
!!provisoire pour changement interface!!!      ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,1)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRIDRYH(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4))
!!provisoire pour changement interface!!!      ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIDRYH(:) / ZRIT(:)        ! QIDRYH
!!provisoire pour changement interface!!!      ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,2)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRSDRYH(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5))
!!provisoire pour changement interface!!!      ZWQ5(:,3) = XCOEF_RQ_S * ZQST(:) * ZRSDRYH(:) / ZRST(:)        ! QSDRYH
!!provisoire pour changement interface!!!      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,3)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRGDRYH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6))
!!provisoire pour changement interface!!!      ZWQ5(:,5) = XCOEF_RQ_G * ZQGT(:) * ZRGDRYH(:) / ZRGT(:)        ! QGDRYH
!!provisoire pour changement interface!!!      ZWQ5(:,5) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQGS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQGS(:) = ZQGS(:) - ZWQ5(:,5)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,5)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    WHERE (ZRRDRYH(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3))
!!provisoire pour changement interface!!!      ZWQ5(:,4) = XCOEF_RQ_R * ZQRT(:) * ZRRDRYH(:) / ZRRT(:)        ! QRDRYH
!!provisoire pour changement interface!!!      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
!!provisoire pour changement interface!!!      ZQHS(:) = ZQHS(:) + ZWQ5(:,4)
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      7.3   conversion of hail into graupel (qdryhg)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRDRYHG, ZRHT, ZQHT, PTSTEP, &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(7), XQTMIN(7), XCOEF_RQ_H, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQHS, ZQGS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'DRYH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'DRYH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'DRYH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 )%PTR%ADD_PHY(D, 'DRYH', &
!!provisoire pour changement interface!!!                             Unpack( -zwq5(:, 3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 )%PTR%ADD_PHY(D, 'DRYH', &
!!provisoire pour changement interface!!!                             Unpack( (-zwq5(:, 5) - zwq(:)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 )%PTR%ADD_PHY(D, 'DRYH',                      &
!!provisoire pour changement interface!!!                             Unpack( ( zwq5(:, 1) + zwq5(:, 2) + zwq5(:, 3) + zwq5(:, 4) + zwq5(:, 5) + zwq(:) ) &
!!provisoire pour changement interface!!!                                       * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      7.4   melting of hail (qhmltr)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    CALL COMPUTE_CHARGE_TRANSFER (ZRHMLTR, ZRHT, ZQHT, PTSTEP, &
!!provisoire pour changement interface!!!                                  XRTMIN_ELEC(7), XQTMIN(7), XCOEF_RQ_H, &
!!provisoire pour changement interface!!!                                  ZWQ, ZQHS, ZQRS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'HMLT', &
!!provisoire pour changement interface!!!                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 )%PTR%ADD_PHY(D, 'HMLT', &
!!provisoire pour changement interface!!!                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  END IF  ! end if krr==7
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      8.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_i
!!provisoire pour changement interface!!!!             ---------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      8.1   Bergeron-Findeisen effect (qcberi)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  CALL COMPUTE_CHARGE_TRANSFER (ZRCBERI, ZRCT, ZQCT, PTSTEP, &
!!provisoire pour changement interface!!!                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
!!provisoire pour changement interface!!!                                ZWQ, ZQCS, ZQIS)
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'BERFI', &
!!provisoire pour changement interface!!!                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'BERFI', &
!!provisoire pour changement interface!!!                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!  end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      9.    COMPUTE THE CHARGE TRANSFER ASSOCIATED WITH THE CORRECTION TERM
!!provisoire pour changement interface!!!!             ---------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  IF (HCLOUD == 'LIMA') THEN
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%INIT_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                              Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%INIT_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                              Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    ZWQ1(:) = 0.
!!provisoire pour changement interface!!!    WHERE (ZRCCORR2(:) .NE. 0. .AND. ZRCT(:) .GT. XRTMIN_ELEC(2))
!!provisoire pour changement interface!!!      ZWQ1(:) = XCOEF_RQ_C * ZQCT(:) * ZRCCORR2(:) / ZRCT(:)
!!provisoire pour changement interface!!!      ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ1(:)) ),ZQCS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQCS(:)  = ZQCS(:)  - ZWQ1(:)
!!provisoire pour changement interface!!!      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ1(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ1(:)/ELECD%XECHARGE )      
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZWQ2(:) = 0.
!!provisoire pour changement interface!!!    WHERE (ZRRCORR2(:) .NE. 0. .AND. ZRRT(:) .GT. XRTMIN_ELEC(3))
!!provisoire pour changement interface!!!      ZWQ2(:) = XCOEF_RQ_R * ZQRT(:) * ZRRCORR2(:) / ZRRT(:)
!!provisoire pour changement interface!!!      ZWQ2(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ2(:)) ),ZQRS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQRS(:)  = ZQRS(:)  - ZWQ2(:)
!!provisoire pour changement interface!!!      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ2(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ2(:)/ELECD%XECHARGE )      
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!    !
!!provisoire pour changement interface!!!    ZWQ3(:) = 0.
!!provisoire pour changement interface!!!    WHERE (ZRICORR2(:) .NE. 0. .AND. ZRIT(:) .GT. XRTMIN_ELEC(4))
!!provisoire pour changement interface!!!      ZWQ3(:) = XCOEF_RQ_I * ZQIT(:) * ZRICORR2(:) / ZRIT(:)
!!provisoire pour changement interface!!!      ZWQ3(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ3(:)) ),ZQIS(:) )
!!provisoire pour changement interface!!!      !
!!provisoire pour changement interface!!!      ZQIS(:)  = ZQIS(:)  - ZWQ3(:)
!!provisoire pour changement interface!!!      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ3(:)/ELECD%XECHARGE )
!!provisoire pour changement interface!!!      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ3(:)/ELECD%XECHARGE )      
!!provisoire pour changement interface!!!    END WHERE
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!    if ( BUCONF%LBUDGET_SV ) then
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG )%PTR%END_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                             Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND )%PTR%END_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                             Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 )%PTR%ADD_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                             Unpack( zwq1(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 )%PTR%ADD_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                             Unpack( zwq2(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )                     
!!provisoire pour changement interface!!!      CALL TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 )%PTR%ADD_PHY(D, 'CORR2', &
!!provisoire pour changement interface!!!                             Unpack( zwq3(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
!!provisoire pour changement interface!!!    end if    
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      X.    COMPUTE THE SEDIMENTATION SOURCE FOR Q_x
!!provisoire pour changement interface!!!!             ----------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!! The sedimentation for electric charges is computed directly 
!!provisoire pour changement interface!!!! in the microphysics scheme
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      10.   UPDATE VOLUMETRIC CHARGE CONCENTRATIONS
!!provisoire pour changement interface!!!!             ---------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!  DO JL = 1, KMICRO
!!provisoire pour changement interface!!!    PQPIS(II1(JL),II2(JL),II3(JL)) = ZQPIS(JL)
!!provisoire pour changement interface!!!    PQNIS(II1(JL),II2(JL),II3(JL)) = ZQNIS(JL)
!!provisoire pour changement interface!!!    PQCS (II1(JL),II2(JL),II3(JL)) = ZQCS(JL)
!!provisoire pour changement interface!!!    PQRS (II1(JL),II2(JL),II3(JL)) = ZQRS(JL)
!!provisoire pour changement interface!!!    PQIS (II1(JL),II2(JL),II3(JL)) = ZQIS(JL)
!!provisoire pour changement interface!!!    PQSS (II1(JL),II2(JL),II3(JL)) = ZQSS(JL)
!!provisoire pour changement interface!!!    PQGS (II1(JL),II2(JL),II3(JL)) = ZQGS(JL)
!!provisoire pour changement interface!!!  END DO
!!provisoire pour changement interface!!!  IF ( KRR == 7 ) THEN
!!provisoire pour changement interface!!!    DO JL = 1, KMICRO
!!provisoire pour changement interface!!!      PQHS(II1(JL),II2(JL),II3(JL)) = ZQHS(JL)
!!provisoire pour changement interface!!!    END DO
!!provisoire pour changement interface!!!  END IF
!!provisoire pour changement interface!!!END IF  ! end if kmicro>0
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!------------------------------------------------------------------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!!*      11.   DEALLOCATE
!!provisoire pour changement interface!!!!             ----------
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!IF (ALLOCATED( ZDELTALWC )) DEALLOCATE( ZDELTALWC )   
!!provisoire pour changement interface!!!IF (ALLOCATED( ZFT ))       DEALLOCATE( ZFT )
!!provisoire pour changement interface!!!!  
!!provisoire pour changement interface!!!IF (ALLOCATED( ZEW ))          DEALLOCATE( ZEW )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSK ))      DEALLOCATE( ZSAUNSK )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIM ))      DEALLOCATE( ZSAUNIM )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIN ))      DEALLOCATE( ZSAUNIN )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSM ))      DEALLOCATE( ZSAUNSM )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSN ))      DEALLOCATE( ZSAUNSN )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZFQIAGGS ))     DEALLOCATE( ZFQIAGGS )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZFQIDRYGBS ))   DEALLOCATE( ZFQIDRYGBS )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZLBQSDRYGB1S )) DEALLOCATE( ZLBQSDRYGB1S )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZLBQSDRYGB2S )) DEALLOCATE( ZLBQSDRYGB2S )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZLBQSDRYGB3S )) DEALLOCATE( ZLBQSDRYGB3S )
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!IF (ALLOCATED( ZDQ ))        DEALLOCATE( ZDQ )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZRAR ))       DEALLOCATE( ZRAR )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZDQ_IS ))     DEALLOCATE( ZDQ_IS )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIM_IS )) DEALLOCATE( ZSAUNIM_IS )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIN_IS )) DEALLOCATE( ZSAUNIN_IS )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZDQ_IG ))     DEALLOCATE( ZDQ_IG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIM_IG )) DEALLOCATE( ZSAUNIM_IG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNIN_IG )) DEALLOCATE( ZSAUNIN_IG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZDQ_SG ))     DEALLOCATE( ZDQ_SG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSK_SG )) DEALLOCATE( ZSAUNSK_SG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSM_SG )) DEALLOCATE( ZSAUNSM_SG )
!!provisoire pour changement interface!!!IF (ALLOCATED( ZSAUNSN_SG )) DEALLOCATE( ZSAUNSN_SG )
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!IF (ALLOCATED( ZEFIELDW ))   DEALLOCATE( ZEFIELDW )
!!provisoire pour changement interface!!!!
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRCMLTSR)) DEALLOCATE(ZRCMLTSR)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRICFRR))  DEALLOCATE(ZRICFRR)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRVHENC))  DEALLOCATE(ZRVHENC)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRCHINC))  DEALLOCATE(ZRCHINC)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRVHONH))  DEALLOCATE(ZRVHONH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRRCVRC))  DEALLOCATE(ZRRCVRC)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRICNVI))  DEALLOCATE(ZRICNVI)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRVDEPI))  DEALLOCATE(ZRVDEPI)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRSHMSI))  DEALLOCATE(ZRSHMSI)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRGHMGI))  DEALLOCATE(ZRGHMGI)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRICIBU))  DEALLOCATE(ZRICIBU)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRIRDSF))  DEALLOCATE(ZRIRDSF)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRCCORR2)) DEALLOCATE(ZRCCORR2)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRRCORR2)) DEALLOCATE(ZRRCORR2)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRICORR2)) DEALLOCATE(ZRICORR2)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRWETGH))  DEALLOCATE(ZRWETGH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRCWETH))  DEALLOCATE(ZRCWETH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRIWETH))  DEALLOCATE(ZRIWETH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRSWETH))  DEALLOCATE(ZRSWETH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRGWETH))  DEALLOCATE(ZRGWETH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRRWETH))  DEALLOCATE(ZRRWETH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRCDRYH))  DEALLOCATE(ZRCDRYH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRRDRYH))  DEALLOCATE(ZRRDRYH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRIDRYH))  DEALLOCATE(ZRIDRYH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRSDRYH))  DEALLOCATE(ZRSDRYH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRGDRYH))  DEALLOCATE(ZRGDRYH)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRHMLTR))  DEALLOCATE(ZRHMLTR)
!!provisoire pour changement interface!!!IF (ALLOCATED(ZRDRYHG))  DEALLOCATE(ZRDRYHG)
!
!------------------------------------------------------------------
END ASSOCIATE
!
CONTAINS
!
! - routines to initialize the non-inductive charging
! - routines to compute the non-inductive charging
! - various useful routines
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_GARDI()
! ##################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Gardiner et al. (1985)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
!
!*     1.      COMPUTE f(DeltaT) AND (LWC - LWC_crit)
!              --------------------------------------
!
GELEC(:,:) = .FALSE.
!
ZDELTALWC(:) = 0.
ZFT(:) = 0.
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) < CST%XTT
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
WHERE (GELEC(:,4))
  ! f(DeltaT)
  ZFT(:) = - 1.7E-5 * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))**3   &
           - 0.003  * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))**2   &
           - 0.05   * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))      &
           + 0.13
  !
  ! LWC - LWC_crit
  ZDELTALWC(:) = (ZRCT(:) * ZRHODREF(:) * 1.E3) - XLWCC   ! (g m^-3)
ENDWHERE
!
END SUBROUTINE ELEC_INIT_NOIND_GARDI
!
!-----------------------------------------------------------------
!
! ################################
  SUBROUTINE ELEC_INIT_NOIND_EWC()
! ################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Saunders et al. (1991)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
!
!*      1.      PARAMETERS FOR POSITIVE NI CHARGING
!               -----------------------------------
!
GELEC(:,:) = .FALSE.
ZDQ(:) = 0.
ZEW(:) = 0.
!
! positive case is the default value
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
  ZFQIAGGS(:)   = XFQIAGGSP
  ZFQIDRYGBS(:) = XFQIDRYGBSP
ELSE IF (CNI_CHARGING == 'TEEWC') THEN
  ZFQIAGGS(:)   = XFQIAGGSP_TAK
  ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
END IF
ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
ZSAUNIM(:) = XIMP     !3.76
ZSAUNIN(:) = XINP     !2.5
ZSAUNSK(:) = XSKP     !52.8
ZSAUNSM(:) = XSMP     !0.44
ZSAUNSN(:) = XSNP     !2.5
!
!
!*      2.      PARAMETERS FOR NEGATIVE NI CHARGING
!               -----------------------------------
!
! Mansell et al. (2005, JGR): droplet collection efficiency of the graupel ~ 0.6-1.0
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZEW(:) = 0.8 * ZRCT(:) * ZRHODREF(:) * 1.E3   ! (g m^-3)
END WHERE
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZEW, ZZT, XSAUNDER, ZDQ)
    !
    WHERE (ZDQ(:) < 0.)
      ZFQIAGGS(:)   = XFQIAGGSN
      ZFQIDRYGBS(:) = XFQIDRYGBSN
    END WHERE
  ELSE IF (CNI_CHARGING == 'TEEWC') THEN
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZEW, ZZT, XTAKA_TM, ZDQ)
    !
    WHERE (ZDQ(:) < 0.)
      ZFQIAGGS(:)   = XFQIAGGSN_TAK
      ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
    END WHERE
  END IF
!
! value of the parameters for the negative case
  WHERE (ZDQ(:) < 0.)
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ZSAUNIM(:) = XIMN      !2.54
    ZSAUNIN(:) = XINN      !2.8
    ZSAUNSK(:) = XSKN      !24.
    ZSAUNSM(:) = XSMN      !0.5
    ZSAUNSN(:) = XSNN      !2.8
  ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_EWC
!
!------------------------------------------------------------------
!
! ################################
  SUBROUTINE ELEC_INIT_NOIND_RAR()
! ################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Saunders and Peck (1998)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
! local variables
REAL, DIMENSION(KMICRO) :: ZRAR_CRIT     ! critical rime accretion rate
REAL, DIMENSION(KMICRO) :: ZVGMEAN, &    ! mean velocity of graupel 
                           ZVSMEAN       ! mean velocity of snow
!
!
!*      1.      COMPUTE THE CRITICAL RIME ACCRETION RATE
!               ----------------------------------------
!
ZRAR_CRIT(:) = 0.
!
IF (CNI_CHARGING == 'SAP98') THEN
!
  WHERE (ZZT(:) <= CST%XTT .AND. ZZT(:) >= (CST%XTT - 23.7)) ! Original from SAP98
    ZRAR_CRIT(:) = 1.0 + 7.93E-2 * (ZZT(:) - CST%XTT) +    &
                         4.48E-2 * (ZZT(:) - CST%XTT)**2 + &
                         7.48E-3 * (ZZT(:) - CST%XTT)**3 + &
                         5.47E-4 * (ZZT(:) - CST%XTT)**4 + &
                         1.67E-5 * (ZZT(:) - CST%XTT)**5 + &
                         1.76E-7 * (ZZT(:) - CST%XTT)**6
  END WHERE
  !
  WHERE (ZZT(:) < (CST%XTT - 23.7) .AND. ZZT(:) > (CST%XTT - 40.)) ! Added by Mansell
    ZRAR_CRIT(:) = 3.4 * (1.0 - (ABS(ZZT(:) - CST%XTT + 23.7) / & ! et al. (2005)
                   (-23.7 + 40.))**3.)
  END WHERE
  !
  GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT
!
ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
!
  WHERE (ZZT(:) > (CST%XTT - 10.7))
    ZRAR_CRIT(:) = 0.66
  END WHERE
  WHERE (ZZT(:) <= (CST%XTT - 10.7) .AND. ZZT(:) >= (CST%XTT - 23.7))
    ZRAR_CRIT(:) = -1.47 - 0.2 * (ZZT(:) - CST%XTT)
  END WHERE
  WHERE (ZZT(:) < (CST%XTT - 23.7) .AND. ZZT(:) > (CST%XTT - 40.))
    ZRAR_CRIT(:) = 3.3
  END WHERE
  !
  GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
               ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
!
ELSE IF (CNI_CHARGING == 'TERAR') THEN
!
  GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT
END IF
!
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
!
!*      2.      INITIALIZATION FOR ICE CRYSTAL - GRAUPEL COLLISIONS
!               ---------------------------------------------------
!
ZDQ_IG(:) = 0.
!
! positive case is the default value
ZSAUNIM_IG(:) = XIMP
ZSAUNIN_IG(:) = XINP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
ZVGMEAN(:) = 0.
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZVGMEAN(:) = XVGCOEF * ZRHOCOR(:) * ZLBDAG(:)**(-ZDG)
  ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVGMEAN(:) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,2)
!
IF (COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for I-G collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_IG)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_IG)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_IG)
    !
  END IF
    !
  WHERE (ZDQ_IG(:) < 0.)
    ZSAUNIM_IG(:) = XIMN
    ZSAUNIN_IG(:) = XINN
  ENDWHERE
ENDIF
!
!
!*      3.      INITIALIZATION FOR ICE CRYSTAL - SNOW COLLISIONS
!               ------------------------------------------------
!
ZDQ_IS(:)  = 0.
!
! positive case is the default value
ZSAUNIM_IS(:) = XIMP
ZSAUNIN_IS(:) = XINP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
ZVSMEAN(:) = 0.
WHERE (ZLBDAS(:) > 0. .AND. ZRCT(:) > 0.)
  ZVSMEAN(:) = XVSCOEF * ZRHOCOR(:) * ZLBDAS(:)**(-ZDS)
  ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVSMEAN(:) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,1)
!
IF (COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for I-S collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_IS)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_IS)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_IS)
    !
  END IF
  !
  WHERE (ZDQ_IS(:) < 0.)
    ZSAUNIM_IS(:) = XIMN
    ZSAUNIN_IS(:) = XINN
  ENDWHERE
ENDIF
!
!
!*      4.      INITIALIZATION FOR GRAUPEL - SNOW COLLISIONS
!               --------------------------------------------
!
ZDQ_SG(:)  = 0.
!
! positive case is the default value
ZSAUNSK_SG(:) = XSKP
ZSAUNSM_SG(:) = XSMP
ZSAUNSN_SG(:) = XSNP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
WHERE (ZVSMEAN(:) > 0. .AND. ZVGMEAN(:) > 0.)
  ZRAR(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * ABS(ZVGMEAN(:) - ZVSMEAN(:)) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,3)
!
IF( COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for S-G collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_SG)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_SG)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_SG)
    !
  END IF
  !
  WHERE (ZDQ_SG(:) < 0.)
    ZSAUNSK_SG(:) = XSKN
    ZSAUNSM_SG(:) = XSMN
    ZSAUNSN_SG(:) = XSNN
  ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_RAR
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_TAKAH()
! ################################## 
!
! Purpose : initialization for the non-inductive charging process
!           following Takahashi (1978)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE f(T, LWC)
!               -----------------
!
ZDQ(:) = 0.
!
ZEW(:) = ZRCT(:) * ZRHODREF(:) * 1.E3      ! (g m^-3)
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                        ZEW, ZZT, XMANSELL, ZDQ)        
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_TAKAH
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_TEEWC()
! ##################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Tsenova and Mitzeva (2009)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
!
!*      1.      PARAMETERS FOR POSITIVE NI CHARGING
!               -----------------------------------
!
GELEC(:,:) = .FALSE.
ZDQ(:) = 0.
!
! positive case is the default value
ZFQIAGGS(:)   = XFQIAGGSP_TAK
ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
ZSAUNIM(:) = XIMP     !3.76
ZSAUNIN(:) = XINP     !2.5
ZSAUNSK(:) = XSKP_TAK !6.5
ZSAUNSM(:) = XSMP     !0.44
ZSAUNSN(:) = XSNP     !2.5
!
!
!*      2.      PARAMETERS FOR NEGATIVE NI CHARGING
!               -----------------------------------   
!
! Compute the effective water content
ZEW(:) = 0.
!
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZEW(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * 1.E3 
END WHERE
!
GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND.  &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                        ZEW, ZZT, XTAKA_TM, ZDQ)        
!
  WHERE (ZDQ(:) < 0.)
    ZFQIAGGS(:)   = XFQIAGGSN_TAK
    ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ZSAUNIM(:) = XIMN      !2.54
    ZSAUNIN(:) = XINN      !2.8
    ZSAUNSK(:) = XSKN_TAK  !2.0
    ZSAUNSM(:) = XSMN      !0.5
    ZSAUNSN(:) = XSNN      !2.8
 ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_TEEWC
!
!------------------------------------------------------------------
!
! #################################################################
  SUBROUTINE ELEC_INI_NI_SAP98(KMICRO, OMASK, PRAR, PRAR_CRIT, PDQ)
! #################################################################
!
!
! Purpose : compute dQ(RAR,T) in the parameterization of Saunders and Peck (1998)
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
INTEGER,                    INTENT(IN)    :: KMICRO
LOGICAL, DIMENSION(KMICRO), INTENT(IN)    :: OMASK
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PRAR       ! Rime accretion rate 
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PRAR_CRIT  ! Critical rime accretion rate
REAL,    DIMENSION(KMICRO), INTENT(INOUT) :: PDQ        ! interpolated dQ
!
!
!*      1.      COMPUTE dQ(RAR, T)
!    ------------------
!
PDQ(:) = 0.
!
! positive region : Mansell et al., 2005
WHERE (OMASK(:) .AND. PRAR(:) > PRAR_CRIT(:))
  PDQ(:) = MAX(0., 6.74 * (PRAR(:) - PRAR_CRIT(:)) * 1.E-15)
ENDWHERE
!
! negative region : Mansell et al. 2005
WHERE (OMASK(:) .AND. PRAR(:) < PRAR_CRIT(:))
  PDQ(:) = MIN(0., 3.9 * (PRAR_CRIT(:) - 0.1) *                   &
                  (4.0 * ((PRAR(:) - (PRAR_CRIT(:) + 0.1) / 2.) / &
                         (PRAR_CRIT(:) - 0.1))**2 - 1.) * 1.E-15)
ENDWHERE
!
END SUBROUTINE ELEC_INI_NI_SAP98
!
!------------------------------------------------------------------
!
! #################################################################
  SUBROUTINE INTERP_DQ_TABLE (KMICRO, KIND_TEMP, KIND_LWC, OMASK, &
                              PLIQ, PTEMP, PTABLE, PDQ)
! ################################################################# 
!
!
! Purpose : interpolate dQ from a lookup table at each gridpoint
!
!
!*  0.  DECLARATIONS
!    ------------
!
IMPLICIT NONE
!
INTEGER,                    INTENT(IN)    :: KMICRO
INTEGER,                    INTENT(IN)    :: KIND_TEMP, KIND_LWC
LOGICAL, DIMENSION(KMICRO), INTENT(IN)    :: OMASK
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PLIQ   ! effective water content or rime accretion rate
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PTEMP  ! temperature
REAL,    DIMENSION(KIND_LWC+1,KIND_TEMP+1), INTENT(IN)    :: PTABLE  ! lookup table
REAL,    DIMENSION(KMICRO), INTENT(INOUT) :: PDQ    ! interpolated dQ
!
! declaration of local variables
INTEGER                 :: IGAUX
REAL,    DIMENSION(:), ALLOCATABLE :: ZDQ_INTERP
REAL,    DIMENSION(:), ALLOCATABLE :: ZVECT1, ZVECT2
INTEGER, DIMENSION(:), ALLOCATABLE :: IVECT1, IVECT2
!
!
!*  1.      FIND THE INDEXES FOR RAR/EW AND T
!    ---------------------------------
!
PDQ(:) = 0.
!
IGAUX = 0
DO II = 1, KMICRO
  IF (OMASK(II)) THEN
    IGAUX = IGAUX + 1
    I1(IGAUX) = II
  END IF 
END DO
!
IF (IGAUX > 0) THEN
  ALLOCATE(ZDQ_INTERP(IGAUX))
  ALLOCATE(ZVECT1(IGAUX))
  ALLOCATE(ZVECT2(IGAUX))
  ALLOCATE(IVECT1(IGAUX))
  ALLOCATE(IVECT2(IGAUX))
  ZDQ_INTERP(:) = 0.
  IVECT1(:) = 0
  IVECT2(:) = 0
!
  DO II = 1, IGAUX
    ZVECT1(II) = PTEMP(I1(II))
    ZVECT2(II) = PLIQ(I1(II))
    ZDQ_INTERP(II) = PDQ(I1(II))
  END DO
!
! Temperature index (0C --> -40C)
  ZVECT1(1:IGAUX) = MAX( 1.00001, MIN( REAL(KIND_TEMP)-0.00001, &
                                    (ZVECT1(1:IGAUX) - CST%XTT - 1.)/(-1.) ) )
  IVECT1(1:IGAUX) = INT( ZVECT1(1:IGAUX) )
  ZVECT1(1:IGAUX) = ZVECT1(1:IGAUX) - REAL(IVECT1(1:IGAUX))
!
! LWC index (0.01 g.m^-3 --> 10 g.m^-3)
  WHERE (ZVECT2(:) >= 0.01 .AND. ZVECT2(:) < 0.1)
    ZVECT2(:) = MAX( 1.00001, MIN( REAL(10)-0.00001, &
                                  ZVECT2(:) * 100. ))
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
  WHERE (ZVECT2(:) >= 0.1 .AND. ZVECT2(:) < 1. .AND. IVECT2(:) == 0)
    ZVECT2(:) = MAX( 10.00001, MIN( REAL(19)-0.00001, &
                                   ZVECT2(:) * 10. + 9. ) )
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
  WHERE ((ZVECT2(:) >= 1.) .AND. ZVECT2(:) <= 10.)
    ZVECT2(:) = MAX( 19.00001, MIN( REAL(KIND_LWC)-0.00001, &
                                   ZVECT2(:) + 18. ) )
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
!
!*  2.      INTERPOLATE dQ(RAR or EW,T)
!    ---------------------------
!
  ZDQ_INTERP(:) = BI_LIN_INTP_V( PTABLE, IVECT2, IVECT1, ZVECT2, ZVECT1, &
                                 IGAUX )
!
  DO II = 1, IGAUX
    PDQ(I1(II)) = ZDQ_INTERP(II)
  END DO
END IF
!
DEALLOCATE(ZDQ_INTERP)
DEALLOCATE(ZVECT1)
DEALLOCATE(ZVECT2)
DEALLOCATE(IVECT1)
DEALLOCATE(IVECT2)
!
END SUBROUTINE INTERP_DQ_TABLE
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_IAGGS_B()
! #########################
!
!
! Purpose : compute charge separation process during the collision
!           between ice and snow
!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!
!*      1.     COMPUTE THE COLLISION EFFICIENCY
!              --------------------------------
!
ZQCOLIS(:) = ZCOLIS * EXP(ZCOLEXIS * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.     COMPUTE THE RATE OF SEPARATED CHARGE
!              ------------------------------------
!
!*      2.1    Charging process following Helsdon and Farley (1987)
!
IF (CNI_CHARGING == 'HELFA') THEN
  !
  WHERE (ZCIT(:) > 0.0 .AND. &
         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5))
    ZWQ_NI(:) = XFQIAGGSBH * ZRIAGGS(:) * ZCIT(:) / ZRIT(:)
    ZWQ_NI(:) = ZWQ_NI(:) * (1. - ZQCOLIS(:)) / ZQCOLIS(:) 
!
! Temperature dependance of the charge transferred
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  END WHERE
!
ELSE
!
!
!*      2.2    Charging process following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    WHERE (GELEC(:,1) .AND. ZDELTALWC(:) > 0.)
      ZWQ_NI(:) = XFQIAGGSBG  * (1 - ZQCOLIS(:)) *         &
                  ZRHODREF(:)**(-4. * ZCEXVT + 4. / ZBI) * &
                  ZCIT(:)**(1 - 4. / ZBI) *                &
                  ZDELTALWC(:) * ZFT(:) *                  &
                  ZCST(:) * ZLBDAS(:)**(-2. - 4. * ZDS) *  &
                  (ZAI * MOMG(ZALPHAI, ZNUI, ZBI) /        &
                  ZRIT(:))**(-4 / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.3    Charging process based on EW: SAUN1/SAUN2, TEEWC
!*             following Saunders et al. (1991), Takahashi via Tsenova and Mitzeva (2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,1) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBS * (1 - ZQCOLIS(:)) *                       &
                  ZRHOCOR(:)**(1 + ZSAUNIN(:)) *                        &
                  ZFQIAGGS(:) * ZDQ(:) *                                &
                  ZCIT(:)**(1 - ZSAUNIM(:) / ZBI) *                     &
                  ZCST(:) * ZLBDAS(:)**(-2.- ZDS * (1. + ZSAUNIN(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.4    Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*             following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies) 
!*             or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    IF (CNI_CHARGING /= 'TERAR') THEN
      ZFQIAGGS(:) = XFQIAGGSP
      WHERE (ZDQ_IS(:) < 0.)
        ZFQIAGGS(:) = XFQIAGGSN
      ENDWHERE
    ELSE
      ZFQIAGGS(:) = XFQIAGGSP_TAK
      WHERE (ZDQ_IS(:) < 0.)
        ZFQIAGGS(:) = XFQIAGGSN_TAK
      ENDWHERE
    ENDIF
!
    WHERE (GELEC(:,1) .AND. ZDQ_IS(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBS * (1 - ZQCOLIS(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IS(:)) *                        &
                  ZFQIAGGS(:) * ZDQ_IS(:) *                                &
                  ZCIT(:)**(1 - ZSAUNIM_IS(:) / ZBI) *                     &
                  ZCST(:) * ZLBDAS(:)**(-2.- ZDS * (1. + ZSAUNIN_IS(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM_IS(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.5    Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    WHERE (GELEC(:,1) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBT1 * (1.0 - ZQCOLIS(:)) * ZRHOCOR(:) *         &
                  ZCIT(:) * ZCST(:) * ZDQ(:) *                            &
                  MIN( XFQIAGGSBT2 / (ZLBDAS(:)**(2. + ZDS)) ,            &
                       XFQIAGGSBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./ZBI) * &
                       ZRIT(:)**(2. / ZBI) /                              &
                      (ZCIT(:)**(2. / ZBI) * ZLBDAS(:)**(2. + 2. * ZDS)))
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_IS
  WHERE (ZWQ_NI(:) .NE. 0.)
    ZLIMIT(:) = XLIM_NI_IS * ZRIAGGS(:) * ZCIT(:) * &
                (1 - ZQCOLIS(:)) / (ZRIT(:) * ZQCOLIS(:))
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
ZQSS(:) = ZQSS(:) + ZWQ_NI(:)
ZQIS(:) = ZQIS(:) - ZWQ_NI(:)
!
END SUBROUTINE ELEC_IAGGS_B
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_IDRYG_B()
! #########################
!
!
!   Purpose : compute charge separation process during the dry collision
!             between ice and graupeln
!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!
!*      1.     COMPUTE THE COLLISION EFFICIENCY
!              --------------------------------
!
ZQCOLIG(:) = ZCOLIG * EXP(ZCOLEXIG * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.     COMPUTE THE RATE OF SEPARATED CHARGE
!              ------------------------------------
!
!*      2.1    Charging process following Helsdon and Farley (1987)
!
IF (CNI_CHARGING == 'HELFA') THEN
  !
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZCIT(:) > 0. .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ_NI(:) = XHIDRYG * ZRIDRYG(:) * ZCIT(:) / ZRIT(:)    
    ZWQ_NI(:) = ZWQ_NI(:) * (1. - ZQCOLIG(:)) / ZQCOLIG(:)        ! QIDRYG_boun
!
! Temperature dependance of the charge transfered
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  END WHERE
!
ELSE
!
!
!*      2.2    Charging process following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    WHERE (GELEC(:,2) .AND. ZDELTALWC(:) > 0.) 
      ZWQ_NI(:) = XFQIDRYGBG * XLBQIDRYGBG * (1 - ZQCOLIG(:)) * &
                  ZRHODREF(:)**(-4. * ZCEXVT + 4. / ZBI) *   &
                  ZCIT(:)**(1 - 4. / ZBI) *                  &
                  ZDELTALWC(:) * ZFT(:) *                    &
                  ZCGT(:) * ZLBDAG(:)**(-2. - 4. * ZDG) *    &
                  (ZAI * MOMG(ZALPHAI, ZNUI, ZBI) /          &
                  ZRIT(:))**(-4 / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.3    Charging process based on EW: SAUN1/SAUN2, TEEWC following
!*             Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,2) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBS * (1. - ZQCOLIG(:)) *                       &
                  ZRHOCOR(:)**(1. + ZSAUNIN(:)) *                        &
                  ZFQIDRYGBS(:) * ZDQ(:) *                               &
                  ZCIT(:)**(1. - ZSAUNIM(:) / ZBI) *                     &
                  ZCGT(:) * ZLBDAG(:)**(-2. - ZDG * (1. + ZSAUNIN(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.4    Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*             following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies)
!*             or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    IF (CNI_CHARGING /= 'TERAR') THEN
       ZFQIDRYGBS(:) = XFQIDRYGBSP
       WHERE (ZDQ_IG(:) < 0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN
       ENDWHERE
    ELSE
       ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
       WHERE (ZDQ_IG(:) <0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
       ENDWHERE
    END IF
!
    WHERE (GELEC(:,2) .AND. ZDQ_IG(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBS * (1. - ZQCOLIG(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IG(:)) *                         &
                  ZFQIDRYGBS(:) * ZDQ_IG(:) *                               &
                  ZCIT(:)**(1 - ZSAUNIM_IG(:) / ZBI) *                      &
                  ZCGT(:) * ZLBDAG(:)**(-2. - ZDG * (1. + ZSAUNIN_IG(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM_IG(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.5     Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    WHERE (GELEC(:,2) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBT1 * (1. - ZQCOLIG(:)) * ZRHOCOR(:) *          &
                  ZCIT(:) * ZCGT(:) * ZDQ(:) *                            &
                  MIN( XFQIDRYGBT2 / (ZLBDAG(:)**(2. + ZDG)),             &
                       XFQIDRYGBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./ZBI) * &
                       ZRIT(:)**(2. / ZBI) / (ZCIT(:)**(2. / ZBI) *       &
                       ZLBDAG(:)**(2. + 2. * ZDG)) )
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_IG
  WHERE (ZWQ_NI(:) .NE. 0. .AND. ZRIT(:) > 0.)
    ZLIMIT(:) = XLIM_NI_IG * ZRIDRYG(:) * ZCIT(:) * (1 - ZQCOLIG(:)) / &
                (ZRIT(:) * ZQCOLIG(:))
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
WHERE (ZRIDRYG(:) > 0.)
  ZQGS(:) = ZQGS(:) + ZWQ_NI(:)
  ZQIS(:) = ZQIS(:) - ZWQ_NI(:)
END WHERE
!
END SUBROUTINE ELEC_IDRYG_B
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_SDRYG_B()
! #########################
!
!
! Purpose : compute the charge separation during the dry collision
!           between snow and graupel
!
!
!*  0.  DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE THE COLLECTION EFFICIENCY
!               ---------------------------------
!
ZQCOLSG(:) = ZCOLSG * EXP (ZCOLEXSG * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.      COMPUTE THE RATE OF SEPARATED CHARGE
!               ------------------------------------
!
!*      2.1     Charge separation following Helsdon and Farley (1987)
! 
IF (CNI_CHARGING == 'HELFA') THEN
!
  WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZLBDAG(:) > 0. .AND. &
         ZRST(:) > XRTMIN_ELEC(5) .AND. ZLBDAS(:) > 0.)
    ZWQ_NI(:) = ZWQ5(:,5) * XFQSDRYGBH * ZRHODREF(:)**(-ZCEXVT) *    &
               (1. - ZQCOLSG(:)) *                                   &
                ZCST(:) * ZCGT(:) *                                  &
               (XLBQSDRYGB4H * ZLBDAS(:)**(-2.) +                    &
                XLBQSDRYGB5H * ZLBDAS(:)**(-1.) * ZLBDAG(:)**(-1.) + &
                XLBQSDRYGB6H * ZLBDAG(:)**(-2.))
!
! Temperature dependance of the charge transfered
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
ELSE
!
!
!*      2.2     Charge separation following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN 
    WHERE (GELEC(:,3) .AND. ZDELTALWC(:) > 0.) 
      ZWQ_NI(:) = XFQSDRYGBG * (1. - ZQCOLSG(:)) *                     &
                  ZRHODREF(:)**(-4. * ZCEXVT) *                        &
                  ZFT(:) * ZDELTALWC(:) *                              &
                  ZCST(:) * ZCGT(:) *                                  &
                 (XLBQSDRYGB4G * ZLBDAS(:)**(-4.) * ZLBDAG(:)**(-2.) + &
                  XLBQSDRYGB5G * ZLBDAS(:)**(-5.) * ZLBDAG(:)**(-1.) + &
                  XLBQSDRYGB6G * ZLBDAS(:)**(-6.)) *                   &
                  ZWQ5(:,5)
    ENDWHERE
  END IF
!
!
!*      2.3     Charging process based on EW: SAUN1/SAUN2, TEEWC following
!*              Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,3) .AND. ZDQ(:) /= 0.) 
!    ZWQ_NI(:) = ZWQ5(:,6) If graupel gains positive charge ZDQ(:) > 0.
!    ZWQ_NI(:) = ZWQ5(:,7) If graupel gains negative charge ZDQ(:) < 0.
      ZWQ_NI(:) = ZWQ5(:,6) * (0.5 + SIGN(0.5,ZDQ(:))) + &
                  ZWQ5(:,7) * (0.5 - SIGN(0.5,ZDQ(:)))
!
      ZWQ_NI(:) = ZWQ_NI(:) * XFQSDRYGBS * (1. - ZQCOLSG(:)) *                  &
                  ZRHOCOR(:)**(1. + ZSAUNSN(:)) * ZSAUNSK(:) * ZDQ(:) *         &
                  ZCST(:) * ZCGT(:) *                                           &
                ( ZLBQSDRYGB1S(:) / (ZLBDAS(:)**ZSAUNSM(:) * ZLBDAG(:)**2) +    &
                  ZLBQSDRYGB2S(:) / (ZLBDAS(:)**( 1.+ZSAUNSM(:)) * ZLBDAG(:)) + &
                  ZLBQSDRYGB3S(:) /  ZLBDAS(:)**(2.+ZSAUNSM(:)) )
    ENDWHERE
  END IF
!
!
!*      2.4     Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*              following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies)
!*              or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
    WHERE (ZDQ_SG(:) < 0.)
      ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
      ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
      ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ENDWHERE
!
    WHERE (GELEC(:,3) .AND. ZDQ_SG(:) /= 0.)
      ZWQ_NI(:) = ZWQ5(:,6) * (0.5+SIGN(0.5,ZDQ_SG(:))) + &
                  ZWQ5(:,7) * (0.5-SIGN(0.5,ZDQ_SG(:)))
!
      ZWQ_NI(:) = ZWQ_NI(:) * XFQSDRYGBS * (1. - ZQCOLSG(:)) *                       &
                  ZRHOCOR(:)**(1. + ZSAUNSN_SG(:)) * ZSAUNSK_SG(:) * ZDQ_SG(:) *     &
                  ZCST(:) * ZCGT(:) *                                                &
                 (ZLBQSDRYGB1S(:) / (ZLBDAS(:)**ZSAUNSM_SG(:)      * ZLBDAG(:)**2) + &
                  ZLBQSDRYGB2S(:) / (ZLBDAS(:)**(1.+ZSAUNSM_SG(:)) * ZLBDAG(:)) +    &
                  ZLBQSDRYGB3S(:) /  ZLBDAS(:)**(2.+ZSAUNSM_SG(:)) )
    ENDWHERE
  END IF
!
!
!*      2.5     Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN 
    WHERE (GELEC(:,3) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQSDRYGBT1 * (1. - ZQCOLSG(:)) * ZRHOCOR(:) *            &
                  ZCGT(:) * ZCST(:) * ZDQ(:) *                              &
                  MIN(10. * (                                               &
                   ABS(XFQSDRYGBT2 / (ZLBDAG(:)**ZDG * ZLBDAS(:)**2.) -     & 
                       XFQSDRYGBT3 / (ZLBDAS(:)**(2. + ZDS))) +             &
                   ABS(XFQSDRYGBT4 / (ZLBDAG(:)**(2. + ZDG)) -              &
                       XFQSDRYGBT5 / (ZLBDAS(:)**ZDS * ZLBDAG(:)**2.)) +    &
                   ABS(XFQSDRYGBT6 / (ZLBDAG(:)**(1. + ZDG) * ZLBDAS(:)) -  &
                       XFQSDRYGBT7 / (ZLBDAS(:)**(1. + ZDS) * ZLBDAG(:)))), &
                   XFQSDRYGBT8 * ZRHOCOR(:) * ZWQ5(:,5) *                   &
                  (XFQSDRYGBT9  / (ZLBDAS(:)**2. * ZLBDAG(:)**2.) +         &
                   XFQSDRYGBT10 / (ZLBDAS(:)**4.) +                         &
                   XFQSDRYGBT11 / (ZLBDAS(:)**3. * ZLBDAG(:))))
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_SG
  WHERE (ZWQ_NI(:) .NE. 0.)
    ZLIMIT(:) = XLIM_NI_SG * ZWQ5(:,4) * XAUX_LIM *  &
                ZRHOCOR(:) * (1. - ZQCOLSG(:)) *     &
                ZCST(:) * ZCGT(:) *                  &
              ( XAUX_LIM1 / ZLBDAS(:)**2           + &
                XAUX_LIM2 /(ZLBDAS(:) * ZLBDAG(:)) + &
                XAUX_LIM3 / ZLBDAG(:)**2 )
!          
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
WHERE (ZRSDRYG(:) > 0.)
  ZQGS(:) = ZQGS(:) + ZWQ_NI(:)
  ZQSS(:) = ZQSS(:) - ZWQ_NI(:)
END WHERE
!
END SUBROUTINE ELEC_SDRYG_B
!
!------------------------------------------------------------------
!
! #########################################################################
  SUBROUTINE COMPUTE_CHARGE_TRANSFER (PR_RATE, PRXT, PQXT, PTSTEP,        &
                                      PRX_THRESH, PQX_THRESH, PCOEF_RQ_X, &
                                      PQ_RATE, PQXS, PQYS                 )
! #########################################################################
!
!  Purpose : compute the charge transfer rate in proportion of the mass transfer rate
!  x --> y
!  q_rate_xy = r_rate_xy * coef_rq_x * qx_t / rx_t
!
!
!*      0.      DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1     Declaration of dummy arguments
!
REAL, INTENT(IN),    DIMENSION(:)  :: PR_RATE     ! Mass exchange rate from x to y
REAL, INTENT(IN),    DIMENSION(:)  :: PRXT        ! Mixing ratio of x at t
REAL, INTENT(IN),    DIMENSION(:)  :: PQXT        ! Electric charge of x at t
REAL, INTENT(IN)                   :: PRX_THRESH  ! Threshold on mixing ratio
REAL, INTENT(IN)                   :: PQX_THRESH  ! Threshold on electric charge
REAL, INTENT(IN)                   :: PCOEF_RQ_X  ! Coefficient for charge exchange
REAL, INTENT(IN)                   :: PTSTEP      ! Time step
REAL, INTENT(INOUT), DIMENSION(:)  :: PQ_RATE     ! Charge exchange rate from x to y
REAL, INTENT(INOUT), DIMENSION(:)  :: PQXS        ! Electric charge of x - source term
REAL, INTENT(INOUT), DIMENSION(:)  :: PQYS        ! Electric charge of y - source term
!
!
!*      0.2     Declaration of local variables
!
!
PQ_RATE(:) = 0.
!
WHERE (PR_RATE(:) > 0. .AND. &
       PRXT(:) > PRX_THRESH .AND. ABS(PQXT(:)) > PQX_THRESH)
! Compute the charge exchanged during the mass tranfer from species x to y
  PQ_RATE(:) = PCOEF_RQ_X * PR_RATE(:) * PQXT(:) / PRXT(:)
! Limit the charge exchanged to the charge available on x at t
  PQ_RATE(:) = SIGN( MIN( ABS(PQXT(:)/PTSTEP),ABS(PQ_RATE(:)) ), PQXT(:)/PTSTEP )
  !
! Update the source terms of x and y
  PQXS(:) = PQXS(:) - PQ_RATE(:)
  PQYS(:) = PQYS(:) + PQ_RATE(:)
END WHERE
!
END SUBROUTINE COMPUTE_CHARGE_TRANSFER
!
!------------------------------------------------------------------
!
! ###########################################################
  FUNCTION BI_LIN_INTP_V(ZT, KI, KJ, PDX, PDY, KN)  RESULT(PY)
! ###########################################################          
!
! Purpose : 
!
!                 |                   |
! ZT(KI(1),KJ(2))-|-------------------|-ZT(KI(2),KJ(2))
!                 |                   | 
!                 |                   |
!              x2-|-------|y(x1,x2)   |
!                 |       |           |
!              PDY|       |           |
!                 |       |           |
!                 |       |           |
!ZT( KI(1),KJ(1))-|-------------------|-ZT(KI(2),KJ(1))
!                 |  PDX  |x1         |
!                 |                   |
!
!*  0.  DECLARATIONS
!            ------------
!
IMPLICIT NONE
!
!*  0.1  Declaration of dummy arguments
!
INTEGER, INTENT(IN)                 :: KN        ! Size of the result vector
INTEGER, INTENT(IN), DIMENSION(KN)  :: KI        ! Tabulated  coordinate
INTEGER, INTENT(IN), DIMENSION(KN)  :: KJ        ! Tabulated  coordinate
REAL,    INTENT(IN), DIMENSION(:,:) :: ZT        ! Tabulated data
REAL,    INTENT(IN), DIMENSION(KN)  :: PDX, PDY  ! 
!
REAL,                DIMENSION(KN)  :: PY         ! Interpolated value
!
!*  0.2  Declaration of local variables
!
INTEGER :: JL        ! Loop index
!
!
!*  1.  INTERPOLATION
!    -------------
!  
DO JL = 1, KN
  PY(JL) = (1.0 - PDX(JL)) * (1.0 - PDY(JL)) * ZT(KI(JL),  KJ(JL))   + &
            PDX(JL)        * (1.0 - PDY(JL)) * ZT(KI(JL)+1,KJ(JL))   + &
            PDX(JL)        * PDY(JL)         * ZT(KI(JL)+1,KJ(JL)+1) + &
           (1.0 - PDX(JL)) * PDY(JL)         * ZT(KI(JL)  ,KJ(JL)+1)
ENDDO
!
END FUNCTION BI_LIN_INTP_V
!
!------------------------------------------------------------------
!
END SUBROUTINE ELEC_TENDENCIES
END MODULE MODE_ELEC_TENDENCIES
