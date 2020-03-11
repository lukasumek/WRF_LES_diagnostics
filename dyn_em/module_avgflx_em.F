!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_avgflx_em

  USE module_bc
  USE module_model_constants
  USE module_wrf_error

CONTAINS

! cLU start
!-------------------------------------------------------------------------------
  subroutine calc_sgs_and_adv(ids, ide, jds, jde, kds, kde,  &
       & ims, ime, jms, jme, kms, kme,           &
       & its, ite, jts, jte, kts, kte,           &
       & config_flags,                           &
       & muu, muv, mut,                          &          
       & msftx, msfty,                           &
       & msfux, msfuy,                           &
       & msfvx, msfvy,                           &
       & rdx, rdy, fnm, fnp,                     &
       & zx, zy, rdz, rdzw,                      & 
       & dnw, dn,                                &
       & TH_PHY, U_PHY, V_PHY,                   &
       & TKE, W,                                 &
       & ZPM, ZPW,                               &
       & xkmh, xkmv, xkhh, xkhv,                 &
       & cf1, cf2, cf3,                          &
       & moist, rho, hfx, ustar,                 &
       & WMASS,                                  &
       & sfcWTH, sfcUW,                          &
       & sgsUU,  sgsVV,  sgsWW,                  &
       & sgsUV,  sgsUW,  sgsVW,                  &
       & sgsUTH, sgsVTH, sgsWTH,                 &
       & THtotaladvX, THtotaladvY, THtotaladvZ,  & 
       & sgsUTHtend, sgsVTHtend, sgsWTHtend) 

    IMPLICIT NONE

    INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
         ims, ime, jms, jme, kms, kme,  &
         its, ite, jts, jte, kts, kte

    TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(IN) ::  &
         TH_PHY,                                 & ! dry theta at mass--levels
         U_PHY, V_PHY,                           &
         TKE,  W,                                & ! staggered velocities
         ZPM, ZPW,                               & 
         zx, zy, rdz, rdzw,                      & 
         xkmh, xkmv, xkhh, xkhv,                 &
         rho                                

    REAL, INTENT(IN) :: cf1, cf2, cf3, rdx, rdy
  
    REAL, DIMENSION( kms:kme) , INTENT(IN) :: fnm, fnp, dnw, dn

    REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(IN) :: muu, muv, mut, &
         msftx, msfty, msfux, msfuy, msfvx, msfvy,  &
         hfx, ustar

    REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_moist ) , INTENT(IN) :: moist

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(INOUT) ::  &
         WMASS, &
         sgsUU,  sgsVV,  sgsWW,   &
         sgsUV,  sgsUW,  sgsVW,   &
         sgsUTH, sgsVTH, sgsWTH,  &
         THtotaladvX, THtotaladvY, THtotaladvZ,  &
         sgsUTHtend, sgsVTHtend, sgsWTHtend

    REAL, DIMENSION(ims:ime, jms:jme) , INTENT(INOUT) ::  &
         sfcWTH, sfcUW

    ! local vars 
    ! some factors for destaggering of W
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) :: upperfac, lowerfac  
    REAL :: tmpfac1, tmpfac2
    ! gradients
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) :: &
               dTHdX, dTHdY, dTHdZ, dUdX, dUdY, dUdZ, &
               dVdX,  dVdY,  dVdZ,  dWdX, dWdY, dWdZ
    ! some local variables for surface fluxes
    REAL, DIMENSION( ims:ime, jms:jme ) :: cpm_sfc, rho_sfc
    REAL, PARAMETER :: cp = 1004.
    INTEGER :: i,j,k,i_end,j_end
    ! cLU end 
    
    cpm_sfc = 0.0
    rho_sfc = 0.0

    CALL wrf_debug ( 20 , ' starting sgs calc...' )

    i_end = min(ite,ide-1)
    j_end = min(jte,jde-1)

    DO j=jts,j_end
       DO k=kts,min(kte,kde-1)
          DO i=its,i_end
             ! do calculation of weighting factors for destaggering of w-level data
             tmpfac1 = ZPM(i,k,j)  -ZPW(i,k,j)
             tmpfac2 = ZPW(i,k+1,j)-ZPM(i,k,j)
             upperfac(i,k,j) = tmpfac1/(tmpfac1+tmpfac2)
             lowerfac(i,k,j) = tmpfac2/(tmpfac1+tmpfac2)

             WMASS(i,k,j)    = W(i,k+1,j)*upperfac(i,k,j) + W(i,k,j)*lowerfac(i,k,j)
          END DO
       END DO
    END DO
    !  set last leve
    WMASS(its:i_end,kte,jts:j_end) =  W(its:i_end,kte,jts:j_end)

    CALL calc_dVARdX(dTHdX, mut, config_flags, TH_PHY, .false., &
                    msftx, msfty, msfux, msfuy, &
                    msfvx, msfvy, rdx, rdy,     &
                    fnm, fnp, cf1, cf2, cf3,    &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdY(dTHdY, mut, config_flags, TH_PHY,.false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdZ(dTHdZ, config_flags, TH_PHY, ZPM, .false., &
                    dn, dnw, rdz, rdzw, fnm, &
                    fnp, cf1, cf2, cf3,           &
                    ids, ide, jds, jde, kds, kde,            &
                    ims, ime, jms, jme, kms, kme,            &
                    its, ite, jts, jte, kts, kte )

    ! for unstaggered U
    CALL calc_dVARdX(dUdX, mut, config_flags, U_PHY, .false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,           &
                    ims, ime, jms, jme, kms, kme,           &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdY(dUdY, mut, config_flags, U_PHY, .false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdZ(dUdZ, config_flags, U_PHY, ZPM, .false., &
                    dn, dnw, rdz, rdzw, fnm, &
                    fnp, cf1, cf2, cf3,           &
                    ids, ide, jds, jde, kds, kde,            &
                    ims, ime, jms, jme, kms, kme,            &
                    its, ite, jts, jte, kts, kte )
    ! for unstaggered V
    CALL calc_dVARdX(dVdX, mut, config_flags, V_PHY, .false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdY(dVdY, mut, config_flags, V_PHY, .false., &
                    msftx, msfty, msfux, msfuy, &
                    msfvx, msfvy, rdx, rdy, &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdZ(dVdZ, config_flags, V_PHY, ZPM, .false.,&
                    dn, dnw, rdz, rdzw, fnm, &
                    fnp, cf1, cf2, cf3,           &
                    ids, ide, jds, jde, kds, kde,            &
                    ims, ime, jms, jme, kms, kme,            &
                    its, ite, jts, jte, kts, kte )
    ! for W: d/dx and d/dy are staggered on w-levels 
    ! diretion
    CALL calc_dVARdX(dWdX, mut, config_flags, W, .false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy,rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    CALL calc_dVARdY(dWdY, mut, config_flags, W, .false., &
                    msftx, msfty, msfux, msfuy, &
                    msfvx, msfvy, rdx, rdy, &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
    ! dWdZ is on mass levels
    CALL calc_dVARdZ(dWdZ, config_flags, W, ZPW, .true.,       &
                     dn, dnw, rdz, rdzw, fnm, &
                     fnp, cf1, cf2, cf3,           &
                     ids, ide, jds, jde, kds, kde,            &
                     ims, ime, jms, jme, kms, kme,            &
                     its, ite, jts, jte, kts, kte )

    ! extrapolate some gradients which do not have values on lowest level
    dTHdZ(its:i_end,kts,jts:j_end) = dTHdZ(its:i_end,1,jts:j_end)*cf1 + &
             dTHdZ(its:i_end,2,jts:j_end)*cf2 + dTHdZ(its:i_end,3,jts:j_end)*cf3 
    dUdZ(its:i_end,kts,jts:j_end)  = dUdZ(its:i_end,1,jts:j_end)*cf1 + &
             dUdZ(its:i_end,2,jts:j_end)*cf2  + dUdZ(its:i_end,3,jts:j_end)*cf3  
    dVdZ(its:i_end,kts,jts:j_end)  = dVdZ(its:i_end,1,jts:j_end)*cf1 + &
             dVdZ(its:i_end,2,jts:j_end)*cf2  + dVdZ(its:i_end,3,jts:j_end)*cf3

    DO j=jts,j_end
       DO k=kts,min(kte,kde)
          DO i=its,i_end  
             ! calculate the full advection term for dry theta at this time step
             ! for the physical space; these are all be on mass-levels 
             THtotaladvX(i,k,j) = -1.0 * U_PHY(i,k,j) * dTHdX(i,k,j)
             THtotaladvY(i,k,j) = -1.0 * V_PHY(i,k,j) * dTHdY(i,k,j)
             ! vertical theta gradients is on w-levels, so destagger that to get 
             ! to mass levels and use W on mass levels too
             THtotaladvZ(i,k,j) = -0.5 * ( (W(i,k,j)  * dTHdZ(i,k,j)) + &
                                           (W(i,k+1,j)* dTHdZ(i,k+1,j)) )

             ! instantaneous SGS variances following Deardorff (1980)
             sgsUU(i,k,j) = -xkmh(i,k,j) * 2.0 * dUdX(i,k,j) + (TKE(i,k,j)*(2.0/3.0)) ! on mass-levels
             sgsVV(i,k,j) = -xkmh(i,k,j) * 2.0 * dVdY(i,k,j) + (TKE(i,k,j)*(2.0/3.0)) ! on mass-levels
             sgsWW(i,k,j) = -xkmv(i,k,j) * 2.0 * dWdZ(i,k,j) + (TKE(i,k,j)*(2.0/3.0)) ! on mass-levels

             ! instantaneous SGS horizontal momentum flux following Deardorff (1980)
             sgsUV(i,k,j) = -xkmh(i,k,j) * ( dUdY(i,k,j) + dVdX(i,k,j) ) ! on mass-levels

             IF (k .EQ. 1) THEN 
               ! calculate a moist c_p 
               cpm_sfc(i,j) = cp * ( 1. + 0.8*moist(i,kts,j,P_QV) )
               rho_sfc(i,j) = rho(i,kts,j) ! get the air density at the surface
               ! caclulate kinematic heat flux from surface sensible heat flux
               sfcWTH(i,j) = hfx(i,j)/( cpm_sfc(i,j)*rho_sfc(i,j))
               ! calcualte kinematic miomemtum flux 
               ! from ustar = (u'w'^2 + v'w'^2)^0.25
               sfcUW(i,j) = ( (ustar(i,j)**4.0) /2.0)**0.5
             ELSE
               ! eddy diffusivity and eddy viscocity from km_opt=2 are on mass-levels
               ! but horizontal gradients of W and vertical gradient of U/V are w-levels 
               ! so we interpolate excahnge coefficients to w-levels
               sgsUW(i,k,j)  = -( 0.5*(xkmv(i,k,j) + xkmv(i,k-1,j)) ) * ( dUdZ(i,k,j) + dWdX(i,k,j) )
               sgsVW(i,k,j)  = -( 0.5*(xkmv(i,k,j) + xkmv(i,k-1,j)) ) * ( dVdZ(i,k,j) + dWdY(i,k,j) )
             END IF
           END DO
        END DO
     END DO
 
     ! now recalculate the SGS trubulent mixing, but seperately in each 
     ! direction in the physical space when using km_opt=2
     sgsUTHtend(its:ite,kts:kte,jts:jte) = 0. ! reset the instantaneous tendency
     CALL horizontal_diffusion_sX(sgsUTHtend, config_flags, TH_PHY, sgsUTH, &
                                 msftx, msfty, msfux, msfuy,            &
                                 msfvx, msfvy, xkhh, rdx, rdy,          &
                                 fnm, fnp, cf1, cf2, cf3,               &
                                 zx, zy, rdz, rdzw, dnw, dn,            &
                                 ids, ide, jds, jde, kds, kde,          &
                                 ims, ime, jms, jme, kms, kme,          &
                                 its, ite, jts, jte, kts, kte)
     sgsVTHtend(its:ite,kts:kte,jts:jte) = 0. ! reset the instantaneous tendency
     CALL horizontal_diffusion_sY(sgsVTHtend, config_flags, TH_PHY, sgsVTH, &
                                 msftx, msfty, msfux, msfuy,            &
                                 msfvx, msfvy, xkhh, rdx, rdy,          &
                                 fnm, fnp, cf1, cf2, cf3,               &
                                 zx, zy, rdz, rdzw, dnw, dn,            &
                                 ids, ide, jds, jde, kds, kde,          &
                                 ims, ime, jms, jme, kms, kme,          &
                                 its, ite, jts, jte, kts, kte)
     sgsWTHtend(its:ite,kts:kte,jts:jte) = 0. ! reset the instantaneous tendency
     CALL vertical_diffusion_sZ(sgsWTHtend, config_flags, TH_PHY, sgsWTH,  &
                                 xkhv, sfcWTH,                          &
                                 rdz, rdzw, fnm, fnp,                   &
                                 ids, ide, jds, jde, kds, kde,          &
                                 ims, ime, jms, jme, kms, kme,          &
                                 its, ite, jts, jte, kts, kte)

     ! now set the surface fluxes to the staggered arrays of the vertical
     ! momentum fluxes
     ! surface heat flux is set in vertical_diffusion_sZ()
     sgsUW(its:i_end,kts,jts:j_end)  = sfcUW(its:i_end,jts:j_end) 
     sgsVW(its:i_end,kts,jts:j_end)  = sfcUW(its:i_end,jts:j_end) 

     ! reset uppermost level of diagnostics to zero (is neglected)
     sgsUU(its:i_end,kte,jts:j_end)  = 0.
     sgsVV(its:i_end,kte,jts:j_end)  = 0.
     sgsWW(its:i_end,kte,jts:j_end)  = 0.
     sgsUW(its:i_end,kte,jts:j_end)  = 0.
     sgsVW(its:i_end,kte,jts:j_end)  = 0.
     sgsUV(its:i_end,kte,jts:j_end)  = 0.
     sgsUTH(its:i_end,kte,jts:j_end) = 0.
     sgsVTH(its:i_end,kte,jts:j_end) = 0.
     sgsWTH(its:i_end,kte,jts:j_end) = 0.
     THtotaladvX(its:i_end,kte,jts:j_end)  = 0.
     THtotaladvY(its:i_end,kte,jts:j_end)  = 0.
     THtotaladvZ(its:i_end,kte,jts:j_end)  = 0.
     
  end subroutine calc_sgs_and_adv


  subroutine calc_meanadv(ids, ide, jds, jde, kds, kde,           &
       & ims, ime, jms, jme, kms, kme,           &
       & its, ite, jts, jte, kts, kte,           &
       & config_flags,                           &
       & muu, muv, mut,                          &
       & msftx, msfty,                           &
       & msfux, msfuy,                           &
       & msfvx, msfvy,                           &
       & rdx, rdy, fnm, fnp,                     &
       & zx, zy, rdz, rdzw,                      &
       & dn, dnw,                                &
       & cf1, cf2, cf3,                          & 
       & THMEAN, UMEAN, VMEAN, WMEAN, ZPM,       &
       & THmeanadvX, THmeanadvY, THmeanadvZ       )

    IMPLICIT NONE

    INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
         ims, ime, jms, jme, kms, kme,  &
         its, ite, jts, jte, kts, kte

    TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(IN) ::  &
         THMEAN, UMEAN, VMEAN, WMEAN,            &
         ZPM,                                    &
         zx, zy, rdz, rdzw

    REAL, INTENT(IN) :: cf1, cf2, cf3, rdx, rdy

    REAL, DIMENSION( kms:kme) , INTENT(IN) :: fnm, fnp, dnw, dn

    REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(IN) :: muu, muv, mut, &
         msftx, msfty, msfux, msfuy, msfvx, msfvy

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(INOUT) ::  &
         THmeanadvX, THmeanadvY, THmeanadvZ

    INTEGER :: i,j,k,i_end,j_end

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: &
         dTHmeandX, dTHmeandY, dTHmeandZ        

    CALL wrf_debug ( 20 , ' starting meanADV calc...' )

    i_end = min(ite,ide-1)
    j_end = min(jte,jde-1)

   ! now calculate teh mean part of the adv
   ! this will only be valid at the end of an averaging interval 
   ! ignore this for now, because we only output it at the end of the
   ! averageing interval (=history_interval)
   CALL wrf_debug(200,'In avg_flx, before calculating mean floww adv')


   CALL calc_dVARdX(dTHmeandX, mut, config_flags, THMEAN, .false., &
                    msftx, msfty, msfux, msfuy, &
                    msfvx, msfvy, rdx, rdy,     &
                    fnm, fnp, cf1, cf2, cf3,    &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
   CALL calc_dVARdY(dTHmeandY, mut, config_flags, THMEAN,.false., &
                    msftx, msfty, msfux, msfuy,   &
                    msfvx, msfvy, rdx, rdy,       &
                    fnm, fnp, cf1, cf2, cf3, &
                    zx, zy, rdz, rdzw, dnw, dn, &
                    ids, ide, jds, jde, kds, kde,          &
                    ims, ime, jms, jme, kms, kme,          &
                    its, ite, jts, jte, kts, kte )
   CALL calc_dVARdZ(dTHmeandZ, config_flags, THMEAN, ZPM, .false., &
                    dn, dnw, rdz, rdzw, fnm, &
                    fnp, cf1, cf2, cf3,           &
                    ids, ide, jds, jde, kds, kde,            &
                    ims, ime, jms, jme, kms, kme,            &
                    its, ite, jts, jte, kts, kte )

   ! extrapolate to lowest w-level
   dTHmeandZ(its:i_end,kts,jts:j_end) = dTHmeandZ(its:i_end,1,jts:j_end)*cf1 + &
                                        dTHmeandZ(its:i_end,2,jts:j_end)*cf2 + &
                                        dTHmeandZ(its:i_end,3,jts:j_end)*cf3
   ! this is calculated on every time step, but only yields valid results at the
   ! end of an averaging interval. at that time output is done, and variables
   ! are written out. otherwise nothing happens with these variables so this
   ! wont cause a problem (other than marginally invrease computing time)
   DO j=jts,j_end
     DO k=kts,kte-1
       DO i=its,i_end
         THmeanadvX(i,k,j) = -1.0 * UMEAN(i,k,j) * dTHmeandX(i,k,j)
         THmeanadvY(i,k,j) = -1.0 * VMEAN(i,k,j) * dTHmeandY(i,k,j)
         ! WMEAN is on mass levels, so destagger dTHdZ gradient  
         THmeanadvZ(i,k,j) = -1.0 * WMEAN(i,k,j) * &
                             0.5*( dTHmeandZ(i,k+1,j)+dTHmeandZ(i,k,j) )
       END DO
     END DO
  END DO

  end subroutine calc_meanadv
  

  subroutine zero_avgflx(avgflx_rum,avgflx_rvm,avgflx_wwm, &
       & ids, ide, jds, jde, kds, kde,           &
       & ims, ime, jms, jme, kms, kme,           &
       & its, ite, jts, jte, kts, kte, do_cu,    &
       & avgflx_cfu1,avgflx_cfd1,avgflx_dfu1,avgflx_efu1,avgflx_dfd1,avgflx_efd1, & ! cLU start
       & tke_tendm_adv,                      &  
       & tke_tendm_hdiff, tke_tendm_vdiff,                 &
       & ru_tendm_pgf, rv_tendm_pgf, rw_tendm_buoy,        & 
       & ru_tendm_cor, rv_tendm_cor, rw_tendm_cor,         & 
       & ru_tendm_curv, rv_tendm_curv, rw_tendm_curv,      & 
       & rt_tendm_rad, rt_tendm_radlw, rt_tendm_radsw,     & 
       & rt_tendm_mp,                                      & 
       & tke_avg_shear, tke_avg_buoy, tke_avg_diss,        & 
       & Umean,  Vmean,  Wmean, THmean,                    & 
       & UVmean, UWmean, VWmean,                           & 
       & UUmean, VVmean, WWmean,                           & 
       & UTHmean, VTHmean, WTHmean,                        &
       & xkmhMEAN, xkmvMEAN,                               &
       & xkhhMEAN, xkhvMEAN,                               &
       & sgsUUmean,  sgsVVmean,  sgsWWmean,                & 
       & sgsUTHmean, sgsVTHmean, sgsWTHmean,               & 
       & sgsUVmean,  sgsUWmean,  sgsVWmean,                &
       & sfcUWmean,  sfcWTHmean, HFXmean,                  & 
       & sgsUTHtendmean, sgsVTHtendmean, sgsWTHtendmean,   & 
       & meanTHtotaladvX, meanTHtotaladvY, meanTHtotaladvZ )  
    ! cLU end
    IMPLICIT NONE

    INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
         ims, ime, jms, jme, kms, kme,  &
         its, ite, jts, jte, kts, kte

    LOGICAL, INTENT(IN) :: do_cu

    REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: & ! cLU start
         tke_tendm_adv,     & 
         tke_tendm_hdiff, tke_tendm_vdiff, & 
         ru_tendm_pgf,  rv_tendm_pgf,   rw_tendm_buoy,      &  
         ru_tendm_cor,  rv_tendm_cor,   rw_tendm_cor,       &   
         ru_tendm_curv, rv_tendm_curv,  rw_tendm_curv,      &
         rt_tendm_rad,  rt_tendm_radlw, rt_tendm_radsw,     &
         rt_tendm_mp,                                       & 
         tke_avg_shear, tke_avg_buoy, tke_avg_diss,         &
         Umean, Vmean, Wmean, THmean,                       &
         UVmean, UWmean, VWmean,                            &
         UUmean, VVmean, WWmean,                            &
         UTHmean, VTHmean, WTHmean,                         &
         xkmhMEAN, xkmvMEAN,                                &
         xkhhMEAN, xkhvMEAN,                                &
         sgsUUmean, sgsVVmean, sgsWWmean,                   &
         sgsUTHmean, sgsVTHmean, sgsWTHmean,                & 
         sgsUVmean,  sgsUWmean, sgsVWmean,                  &   
         sgsUTHtendmean, sgsVTHtendmean, sgsWTHtendmean,    &
         meanTHtotaladvX, meanTHtotaladvY, meanTHtotaladvZ, &
         avgflx_rum,avgflx_rvm,avgflx_wwm  

    REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT) ::  sfcUWmean, sfcWTHmean, HFXmean
    ! cLU end

    REAL,     OPTIONAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::    &
         avgflx_cfu1,avgflx_cfd1,avgflx_dfu1,avgflx_efu1,avgflx_dfd1,avgflx_efd1

    INTEGER :: i,j,k

    DO j=jts,jte
       DO k=kts,kte
          DO i=its,ite
             ! cLU start 
             IF ( k .EQ. kts ) THEN
                sfcUWmean(i,j)  = 0.
                sfcWTHmean(i,j) = 0. 
                HFXmean(i,j) = 0.
             END IF
 
             ! TKE
             tke_tendm_adv(i,k,j)  = 0.
             tke_tendm_hdiff(i,k,j) = 0.
             tke_tendm_vdiff(i,k,j) = 0.
                          
             ! tendency from pressure and buoyany force 
             ru_tendm_pgf(i,k,j)   = 0.
             rv_tendm_pgf(i,k,j)   = 0. 
             rw_tendm_buoy(i,k,j)  = 0.
             ! tendency from coriolis force 
             ru_tendm_cor(i,k,j)   = 0.
             rv_tendm_cor(i,k,j)   = 0.
             rw_tendm_cor(i,k,j)   = 0.
             ! tendency from curvature effect 
             ru_tendm_curv(i,k,j)  = 0.
             rv_tendm_curv(i,k,j)  = 0.
             rw_tendm_curv(i,k,j)  = 0. 

             ! pot. temperature change from radiation and micro physics
             rt_tendm_rad(i,k,j)   = 0. 
             rt_tendm_radlw(i,k,j) = 0. 
             rt_tendm_radsw(i,k,j) = 0.
             rt_tendm_mp(i,k,j)    = 0.

             ! TKE production and dissipationfrom deardorff (km_opt = 2)
             tke_avg_shear(i,k,j) = 0.
             tke_avg_buoy(i,k,j)  = 0.
             tke_avg_diss(i,k,j)  = 0.

             ! average values U,V,W, dry potential temperature and individual products
             Umean(i,k,j)  = 0.
             Vmean(i,k,j)  = 0.
             Wmean(i,k,j)  = 0.
             THmean(i,k,j) = 0.

             UVmean(i,k,j) = 0. 
             UWmean(i,k,j) = 0.
             VWmean(i,k,j) = 0.
             UUmean(i,k,j) = 0.
             VVmean(i,k,j) = 0.
             WWmean(i,k,j) = 0.
             UTHmean(i,k,j)= 0.
             VTHmean(i,k,j)= 0.
             WTHmean(i,k,j)= 0.

             xkmhMEAN(i,k,j)=0.
             xkmvMEAN(i,k,j)=0.
             xkhhMEAN(i,k,j)=0.
             xkhvMEAN(i,k,j)=0.
 
             ! SGS (co-)variances
             sgsUUmean(i,k,j)  = 0. 
             sgsVVmean(i,k,j)  = 0. 
             sgsWWmean(i,k,j)  = 0.
             sgsUTHmean(i,k,j) = 0.
             sgsVTHmean(i,k,j) = 0.
             sgsWTHmean(i,k,j) = 0.
             sgsUVmean(i,k,j)  = 0. 
             sgsUWmean(i,k,j)  = 0. 
             sgsVWmean(i,k,j)  = 0.   

             ! sgs tendencies from Deardorff (km_op=2)
             sgsUTHtendmean(i,k,j) = 0.
             sgsVTHtendmean(i,k,j) = 0.
             sgsWTHtendmean(i,k,j) = 0.

             ! averaged total advection of dry pot. temperature 
             ! in physical space
             meanTHtotaladvX(i,k,j) = 0.
             meanTHtotaladvY(i,k,j) = 0.
             meanTHtotaladvZ(i,k,j) = 0.  
             ! cLU end

             avgflx_rum(i,k,j) = 0.
             avgflx_rvm(i,k,j) = 0.
             avgflx_wwm(i,k,j) = 0.
          end DO
       end DO
    end DO

    if (do_cu .and. &
         & present(avgflx_cfu1) .and. present(avgflx_cfd1) .and. present(avgflx_dfu1) &
         & .and. present(avgflx_efu1) .and. present(avgflx_dfd1) .and. present(avgflx_efd1) ) then
       DO j=jts,jte
          DO k=kts,kte
             DO i=its,ite
                avgflx_cfu1(i,k,j) = 0.
                avgflx_cfd1(i,k,j) = 0.
                avgflx_dfu1(i,k,j) = 0.
                avgflx_efu1(i,k,j) = 0.
                avgflx_dfd1(i,k,j) = 0.
                avgflx_efd1(i,k,j) = 0.
             end DO
          end DO
       end DO
    end if

    return
  end subroutine zero_avgflx

  subroutine upd_avgflx(avgflx_count,avgflx_rum,avgflx_rvm,avgflx_wwm, &
       & ru_m, rv_m, ww_m, &
       & ids, ide, jds, jde, kds, kde,           &
       & ims, ime, jms, jme, kms, kme,           &
       & its, ite, jts, jte, kts, kte, do_cu,    &
       & cfu1,cfd1,dfu1,efu1,dfd1,efd1,          &
       & avgflx_cfu1,avgflx_cfd1,avgflx_dfu1,avgflx_efu1,avgflx_dfd1,avgflx_efd1, & ! cLU start
       & config_flags,           & 
       & muu, muv, mut,          & ! dry air masses for uncoupling
       & c1, c2, c1f, c2f,       & ! facotrs of vert.coordinate for uncouplint (mass and w-levels)
       & tke_tendm_adv, & ! Mean advection tend
       & tke_tendm_hdiff, tke_tendm_vdiff, &
       & ru_tendm_pgf,  rv_tendm_pgf,   rw_tendm_buoy,   & ! pressure and buoyancy
       & ru_tendm_cor,  rv_tendm_cor,   rw_tendm_cor,    & ! coriolis tendency
       & ru_tendm_curv, rv_tendm_curv,  rw_tendm_curv,   & ! curvature effect
       & rt_tendm_rad,  rt_tendm_radlw, rt_tendm_radsw,  & ! Thetatend from physics
       & rt_tendm_mp,                                    & ! Thetatend from MP
       & tke_avg_shear, tke_avg_buoy, tke_avg_diss,      & ! TKE prod./diss. from deardroff
       & Umean,    Vmean,   Wmean,   THmean,                    &
       & UVmean,   UWmean,  VWmean,                         &
       & UUmean,   VVmean,  WWmean,                         &
       & UTHmean,  VTHmean, WTHmean,                      &
       & xkmhMEAN, xkmvMEAN,                             &
       & xkhhMEAN, xkhvMEAN,                             &
       & sgsUUmean,  sgsVVmean,  sgsWWmean,               & ! SGS variances and TKE
       & sgsUTHmean, sgsVTHmean, sgsWTHmean,                      & ! SGS heat fluxes
       & sgsUVmean,  sgsUWmean,  sgsVWmean,                       & ! SGS mom. fluxes
       & sfcUWmean,  sfcWTHmean, HFXmean,                      & ! SFC fluxes
       & sgsUTHtendmean, sgsVTHtendmean, sgsWTHtendmean,          &
       & tke_tend_adv,                         & ! now instantaneous values: ADV
       & tke_tend_hdiff, tke_tend_vdiff,                 &
       & ru_tend_pgf, rv_tend_pgf, rw_tend_buoy,         & ! pressure and buoyancy
       & ru_tend_cor, rv_tend_cor, rw_tend_cor,          & ! corioils tendency
       & ru_tend_curv, rv_tend_curv, rw_tend_curv,       & ! curvature effect
       & rt_tend_rad, rt_tend_radlw, rt_tend_radsw,      & ! Thetatend from physics
       & rt_tend_mp,                                     & ! Thetatend from MP (diabatic heating)
       & tke_shear, tke_buoy, tke_diss,                  & ! TKE prod./diss. from deardroff
       & U, V, W, WMASS, THETA,        & !
       & Ustag, Vstag,                          &
       & ZPM, ZPW,                              &
       & xkmh, xkmv, xkhh, xkhv,  &
       & cf1, cf2, cf3,           &
       & sgsUU,  sgsVV,  sgsWW,   &
       & sgsUV,  sgsUW,  sgsVW,   &
       & sgsUTH, sgsVTH, sgsWTH,  &
       & sfcUW, sfcWTH, HFX,      &
       & msftx, msfty,                           &
       & msfux, msfuy,                           &
       & msfvx, msfvy,                           &
       & rdx, rdy, fnm, fnp,                     &
       & zx, zy, rdz, rdzw,                      &
       & dnw, dn,                                & 
       & sgsUTHtend, sgsVTHtend, sgsWTHtend,     &
       & THtotaladvX, THtotaladvY, THtotaladvZ,  &
       & meanTHtotaladvX, meanTHtotaladvY, meanTHtotaladvZ )

       ! cLU end 

    IMPLICIT NONE

    INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
         ims, ime, jms, jme, kms, kme,  &
         its, ite, jts, jte, kts, kte

    INTEGER , INTENT(IN)        :: avgflx_count
    LOGICAL, INTENT(IN) :: do_cu
    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(IN) ::  &
         ru_m, rv_m, ww_m,                                      & ! cLU start
         tke_tend_adv,                                &
         tke_tend_hdiff, tke_tend_vdiff,                        &
         ru_tend_pgf,  rv_tend_pgf,  rw_tend_buoy,              &
         ru_tend_cor,  rv_tend_cor,  rw_tend_cor,               &
         ru_tend_curv, rv_tend_curv, rw_tend_curv,              &
         rt_tend_rad, rt_tend_radlw, rt_tend_radsw, rt_tend_mp, &      
         tke_shear, tke_buoy, tke_diss,                         &
         U, V, W, WMASS, THETA,                                 &
         Ustag, Vstag,                    &  
         ZPM, ZPW,                        &
         xkmh, xkmv, xkhh, xkhv,  &
         sgsUU,  sgsVV,  sgsWW,   &
         sgsUV,  sgsUW,  sgsVW,   &
         sgsUTH, sgsVTH, sgsWTH,  &
         zx, zy, rdz, rdzw,       &
         sgsUTHtend, sgsVTHtend, sgsWTHtend, &
         THtotaladvX, THtotaladvY, THtotaladvZ

    REAL , INTENT(IN) ::  cf1, cf2, cf3, rdx, rdy

    TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

    REAL, DIMENSION( kms:kme) , INTENT(IN) :: fnm, fnp, dnw, dn

    REAL, DIMENSION(ims:ime, jms:jme) , INTENT(IN) :: sfcUW, sfcWTH, hfx, & 
         muu, muv, mut,   &
         msftx, msfty, msfux, msfuy, msfvx, msfvy

    ! some factors for uncoupling, taken from phy_prep_part2
    REAL, DIMENSION( kms:kme ) , INTENT(IN) :: c1, c2, c1f, c2f ! uncoupling factors on mass levels and full levels

    REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::    &
         avgflx_rum,avgflx_rvm,avgflx_wwm,                           & ! cLU start
         tke_tendm_adv,                                &
         tke_tendm_hdiff, tke_tendm_vdiff,                           &
         ru_tendm_pgf,  rv_tendm_pgf,   rw_tendm_buoy,               &
         ru_tendm_cor,  rv_tendm_cor,   rw_tendm_cor,                &
         ru_tendm_curv, rv_tendm_curv,  rw_tendm_curv,               &
         rt_tendm_rad,  rt_tendm_radlw, rt_tendm_radsw, rt_tendm_mp, &  
         tke_avg_shear, tke_avg_buoy, tke_avg_diss,                  & 
         Umean, Vmean, Wmean, THmean,                     &
         UVmean, UWmean, VWmean,                          &
         UUmean, VVmean, WWmean,                          &
         UTHmean, VTHmean, WTHmean,                       &
         xkmhMEAN, xkmvMEAN,                              &
         xkhhMEAN, xkhvMEAN,                              &
         sgsUUmean,  sgsVVmean,  sgsWWmean,                  &
         sgsUTHmean, sgsVTHmean, sgsWTHmean,                       &
         sgsUVmean,  sgsUWmean,  sgsVWmean,                  & 
         sgsUTHtendmean, sgsVTHtendmean, sgsWTHtendmean,   &
         meanTHtotaladvX, meanTHtotaladvY, meanTHtotaladvZ

    REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT) :: sfcUWmean, sfcWTHmean, HFXmean
    
    REAL, PARAMETER :: cp = 1004.
    ! cLU end

    REAL,     OPTIONAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN) ::    &
         cfu1,cfd1,dfu1,efu1,dfd1,efd1
    REAL,     OPTIONAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::    &
         avgflx_cfu1,avgflx_cfd1,avgflx_dfu1,avgflx_efu1,avgflx_dfd1,avgflx_efd1

    INTEGER :: i,j,k,i_end,j_end
    REAL :: local_count
   
    i_end = min(ite,ide-1)
    j_end = min(jte,jde-1)
    
    local_count = real(avgflx_count)
    DO j=jts,j_end
       DO k=kts,kte-1
          DO i=its,i_end
             ! cLU start
             IF ( k .EQ. kts) THEN
               ! ustar = (u'w'^2 + v'w'^2)^0.25
               HFXmean(i,j)    = (local_count*HFXmean(i,j)    + hfx(i,j) )   / (local_count+1.)
               sfcWTHmean(i,j) = (local_count*sfcWTHmean(i,j) + sfcWTH(i,j)) / (local_count+1.)
               sfcUWmean(i,j)  = (local_count*sfcUWmean(i,j)  + sfcUW(i,j) ) / (local_count+1.)
             END IF

             ! doing TKE
             tke_tendm_adv(i,k,j)  = (local_count*tke_tendm_adv(i,k,j)  + (tke_tend_adv(i,k,j)  /(c1(k)*mut(i,j)+c2(k)) ) )/(local_count+1.)
             tke_tendm_hdiff(i,k,j)= (local_count*tke_tendm_hdiff(i,k,j)+ (tke_tend_hdiff(i,k,j)/(c1(k)*mut(i,j)+c2(k)) ) )/(local_count+1.)
             tke_tendm_vdiff(i,k,j)= (local_count*tke_tendm_vdiff(i,k,j)+ (tke_tend_vdiff(i,k,j)/(c1(k)*mut(i,j)+c2(k)) ) )/(local_count+1.)

             ! doing pressure gradient and buoyancy
             ru_tendm_pgf(i,k,j) = (local_count*ru_tendm_pgf(i,k,j) + (ru_tend_pgf(i,k,j)/(c1(k)  *muu(i,j)+c2(k) )  ) )/(local_count+1.)
             rv_tendm_pgf(i,k,j) = (local_count*rv_tendm_pgf(i,k,j) + (rv_tend_pgf(i,k,j)/(c1(k)  *muv(i,j)+c2(k) )  ) )/(local_count+1.)
             rw_tendm_buoy(i,k,j)= (local_count*rw_tendm_buoy(i,k,j)+ (rw_tend_buoy(i,k,j)/(c1f(k)*mut(i,j)+c2f(k))  ) )/(local_count+1.)
             ! doing coriolis force 
             ru_tendm_cor(i,k,j) = (local_count*ru_tendm_cor(i,k,j) + (ru_tend_cor(i,k,j)/(c1(k) *muu(i,j)+c2(k))  ) )/(local_count+1.)
             rv_tendm_cor(i,k,j) = (local_count*rv_tendm_cor(i,k,j) + (rv_tend_cor(i,k,j)/(c1(k) *muv(i,j)+c2(k))  ) )/(local_count+1.)
             rw_tendm_cor(i,k,j) = (local_count*rw_tendm_cor(i,k,j) + (rw_tend_cor(i,k,j)/(c1f(k)*mut(i,j)+c2f(k)) ) )/(local_count+1.)
             ! doing curvature effect
             ru_tendm_curv(i,k,j) = (local_count*ru_tendm_curv(i,k,j) + (ru_tend_curv(i,k,j)/(c1(k) *muu(i,j)+c2(k))  ) )/(local_count+1.)
             rv_tendm_curv(i,k,j) = (local_count*rv_tendm_curv(i,k,j) + (rv_tend_curv(i,k,j)/(c1(k) *muv(i,j)+c2(k))  ) )/(local_count+1.)
             rw_tendm_curv(i,k,j) = (local_count*rw_tendm_curv(i,k,j) + (rw_tend_curv(i,k,j)/(c1f(k)*mut(i,j)+c2f(k)) ) )/(local_count+1.)

             ! doing deardroff diss-/prod.
             tke_avg_shear(i,k,j) = (local_count*tke_avg_shear(i,k,j) + (tke_shear(i,k,j) /(c1(k) *mut(i,j)+c2(k)) ) )/(local_count+1.)
             tke_avg_buoy(i,k,j)  = (local_count*tke_avg_buoy(i,k,j)  + (tke_buoy(i,k,j)  /(c1(k) *mut(i,j)+c2(k)) ) )/(local_count+1.)
             tke_avg_diss(i,k,j)  = (local_count*tke_avg_diss(i,k,j)  + (tke_diss(i,k,j)  /(c1(k) *mut(i,j)+c2(k)) ) )/(local_count+1.)

             ! averageing of uncoupled, unstaggered velocity components 
             UMEAN(i,k,j) = (local_count*UMEAN(i,k,j)  + U(i,k,j) )/(local_count+1.)
             VMEAN(i,k,j) = (local_count*VMEAN(i,k,j)  + V(i,k,j) )/(local_count+1.)
             WMEAN(i,k,j) = (local_count*WMEAN(i,k,j)  + WMASS(i,k,j) )/(local_count+1.)
             ! and dry potential temperature             
             THMEAN(i,k,j) = (local_count*THMEAN(i,k,j) + THETA(i,k,j) )/(local_count+1.)

            ! averageing of uncoupled, unstaggered products to derive resolved (co-)variances  
             UVMEAN(i,k,j) = (local_count*UVMEAN(i,k,j) + ( U(i,k,j)*V(i,k,j)     ) )/(local_count+1.)
             UWMEAN(i,k,j) = (local_count*UWMEAN(i,k,j) + ( U(i,k,j)*WMASS(i,k,j) ) )/(local_count+1.)
             VWMEAN(i,k,j) = (local_count*VWMEAN(i,k,j) + ( V(i,k,j)*WMASS(i,k,j) ) )/(local_count+1.)
             UUMEAN(i,k,j) = (local_count*UUMEAN(i,k,j) + ( U(i,k,j)**2           ) )/(local_count+1.)
             VVMEAN(i,k,j) = (local_count*VVMEAN(i,k,j) + ( V(i,k,j)**2           ) )/(local_count+1.)
             WWMEAN(i,k,j) = (local_count*WWMEAN(i,k,j) + ( WMASS(i,k,j)**2       ) )/(local_count+1.)

             UTHMEAN(i,k,j)  = (local_count*UTHMEAN(i,k,j) + ( U(i,k,j)*THETA(i,k,j)     ) )/(local_count+1.)
             VTHMEAN(i,k,j)  = (local_count*VTHMEAN(i,k,j) + ( V(i,k,j)*THETA(i,k,j)     ) )/(local_count+1.)
             WTHMEAN(i,k,j)  = (local_count*WTHMEAN(i,k,j) + ( WMASS(i,k,j)*THETA(i,k,j) ) )/(local_count+1.)

             ! averaging of SGS turbulence exxchange coefficients
             xkmhMEAN(i,k,j) = (local_count*xkmhMEAN(i,k,j) + xkmh(i,k,j) )/(local_count+1.)
             xkmvMEAN(i,k,j) = (local_count*xkmvMEAN(i,k,j) + xkmv(i,k,j) )/(local_count+1.)
             xkhhMEAN(i,k,j) = (local_count*xkhhMEAN(i,k,j) + xkhh(i,k,j) )/(local_count+1.)
             xkhvMEAN(i,k,j) = (local_count*xkhvMEAN(i,k,j) + xkhv(i,k,j) )/(local_count+1.)

             ! doing SGS turbulence 
             sgsUUmean(i,k,j)  = (local_count*sgsUUmean(i,k,j)  + ABS(sgsUU(i,k,j) )) / (local_count+1.)
             sgsVVmean(i,k,j)  = (local_count*sgsVVmean(i,k,j)  + ABS(sgsVV(i,k,j) )) / (local_count+1.)
             sgsWWmean(i,k,j)  = (local_count*sgsWWmean(i,k,j)  + ABS(sgsWW(i,k,j) )) / (local_count+1.)
             sgsUTHmean(i,k,j) = (local_count*sgsUTHmean(i,k,j) + sgsUTH(i,k,j) ) / (local_count+1.)
             sgsVTHmean(i,k,j) = (local_count*sgsVTHmean(i,k,j) + sgsVTH(i,k,j) ) / (local_count+1.)
             sgsWTHmean(i,k,j) = (local_count*sgsWTHmean(i,k,j) + sgsWTH(i,k,j) ) / (local_count+1.)
             sgsUVmean(i,k,j)  = (local_count*sgsUVmean(i,k,j)  + sgsUV(i,k,j) )  / (local_count+1.)
             sgsUWmean(i,k,j)  = (local_count*sgsUWmean(i,k,j)  + sgsUW(i,k,j) )  / (local_count+1.)
             sgsVWmean(i,k,j)  = (local_count*sgsVWmean(i,k,j)  + sgsVW(i,k,j) )  / (local_count+1.)
       
             ! average dry potential temperature tendencies from SGS turblent mixing
             sgsUTHtendmean(i,k,j) = (local_count*sgsUTHtendmean(i,k,j) + sgsUTHtend(i,k,j) ) / (local_count+1.)
             sgsVTHtendmean(i,k,j) = (local_count*sgsVTHtendmean(i,k,j) + sgsVTHtend(i,k,j) ) / (local_count+1.)
             sgsWTHtendmean(i,k,j) = (local_count*sgsWTHtendmean(i,k,j) + sgsWTHtend(i,k,j) ) / (local_count+1.)
 
             ! average tendency coming from full dry potential temperature advection 
             meanTHtotaladvX(i,k,j) = (local_count*meanTHtotaladvX(i,k,j) + THtotaladvX(i,k,j) ) / (local_count+1.)
             meanTHtotaladvY(i,k,j) = (local_count*meanTHtotaladvY(i,k,j) + THtotaladvY(i,k,j) ) / (local_count+1.)
             meanTHtotaladvZ(i,k,j) = (local_count*meanTHtotaladvZ(i,k,j) + THtotaladvZ(i,k,j) ) / (local_count+1.)

             ! averaging uncoupled, unstaggered tendencies from phyiscs parameterizations
             rt_tendm_rad(i,k,j)   = (local_count*rt_tendm_rad(i,k,j)   + rt_tend_rad(i,k,j)   )/(local_count+1.)
             rt_tendm_radlw(i,k,j) = (local_count*rt_tendm_radlw(i,k,j) + rt_tend_radlw(i,k,j) )/(local_count+1.) 
             rt_tendm_radsw(i,k,j) = (local_count*rt_tendm_radsw(i,k,j) + rt_tend_radsw(i,k,j) )/(local_count+1.) 
             rt_tendm_mp(i,k,j)    = (local_count*rt_tendm_mp(i,k,j)    + rt_tend_mp(i,k,j)    )/(local_count+1.)
             ! cLU end
             
             avgflx_rum(i,k,j) = (local_count*avgflx_rum(i,k,j) + ru_m(i,k,j))/(local_count+1.)
             avgflx_rvm(i,k,j) = (local_count*avgflx_rvm(i,k,j) + rv_m(i,k,j))/(local_count+1.)
             avgflx_wwm(i,k,j) = (local_count*avgflx_wwm(i,k,j) + ww_m(i,k,j))/(local_count+1.)
          end DO
       end DO
    end DO

    if (do_cu .and. &
         & present(avgflx_cfu1) .and. present(avgflx_cfd1) .and. present(avgflx_dfu1) &
         & .and. present(avgflx_efu1) .and. present(avgflx_dfd1) .and. present(avgflx_efd1) &
         & .and. present(cfu1) .and. present(cfd1) .and. present(dfu1) &
         & .and. present(efu1) .and. present(dfd1) .and. present(efd1) ) then
       DO j=jts,j_end
          DO k=kts,kte
             DO i=its,i_end
                avgflx_cfu1(i,k,j) = (local_count*avgflx_cfu1(i,k,j) + &
                     & cfu1(i,k,j)) / (local_count+1.)
                avgflx_cfd1(i,k,j) = (local_count*avgflx_cfd1(i,k,j) + &
                     & cfd1(i,k,j)) / (local_count+1.)
                avgflx_dfu1(i,k,j) = (local_count*avgflx_dfu1(i,k,j) + &
                     & dfu1(i,k,j)) / (local_count+1.)
                avgflx_efu1(i,k,j) = (local_count*avgflx_efu1(i,k,j) + &
                     & efu1(i,k,j)) / (local_count+1.)
                avgflx_dfd1(i,k,j) = (local_count*avgflx_dfd1(i,k,j) + &
                     & dfd1(i,k,j)) / (local_count+1.)
                avgflx_efd1(i,k,j) = (local_count*avgflx_efd1(i,k,j) + &
                     & efd1(i,k,j)) / (local_count+1.)
             end DO
          end DO
       end DO
    end if

    return
  end subroutine upd_avgflx


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE calc_dVARdX(dVARdX, mu, config_flags, var, iflux,  &
                       msftx, msfty, msfux, msfuy,               &
                       msfvx, msfvy, rdx, rdy,                   &
                       fnm, fnp, cf1, cf2, cf3,                  &
                       zx, zy, rdz, rdzw, dnw, dn,               &
                       ids, ide, jds, jde, kds, kde,             &
                       ims, ime, jms, jme, kms, kme,             &
                       its, ite, jts, jte, kts, kte)

   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   LOGICAL , INTENT(IN   )  :: iflux  ! .true. : var is a flux
                                      ! .false.: var is no flux

   REAL , INTENT(IN   )           ::        cf1, cf2, cf3

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::     dn
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    dnw

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfux
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfuy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msftx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfty

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   mu

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: dVARdX

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::         &
                                                                     rdz, &
                                                                     rdzw

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::    var, &
                                                                     zx, &
                                                                     zy

   REAL ,                                        INTENT(IN   ) ::    rdx, &
                                                                     rdy

   ! Local data
   INTEGER :: i, j, k, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   REAL , DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1)    ::     H1avg, &
                                                                  H2avg, &
                                                                  H1, &
                                                                  H2
   REAL    :: mrdx, mrdy
   REAL    :: tmpzx, tmpzy, tmpzeta_z, rdzu, rdzv
   INTEGER :: ktes1,ktes2

! End declarations.

!======================================================================

   ktf=MIN(kte,kde-1)

   ktes1=kte-1
   ktes2=kte-2

   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
       config_flags%nested) j_end   = MIN(jde-2,jte)
   IF ( config_flags%periodic_x ) i_start = its
   IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

   ! H1 = partial var over partial x
   DO j = j_start, j_end
   DO k = kts+1, ktf
   DO i = i_start, i_end + 1
      H1avg(i,k,j)=0.5*(fnm(k)*(var(i-1,k  ,j)+var(i,k  ,j))+  &
                        fnp(k)*(var(i-1,k-1,j)+var(i,k-1,j)))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO i = i_start, i_end + 1
      IF (iflux .eqv. .true.) THEN
         H1avg(i,kts  ,j)=0.
         H1avg(i,ktf+1,j)=0.
      ELSEIF (iflux .eqv. .false.) THEN
         H1avg(i,kts  ,j)=0.5*(cf1*var(i  ,1,j)+cf2*var(i  ,2,j)+ &
                               cf3*var(i  ,3,j)+cf1*var(i-1,1,j)+  &
                               cf2*var(i-1,2,j)+cf3*var(i-1,3,j))
      
         H1avg(i,ktf+1,j)=0.5*(var(i,ktes1,j)+(var(i,ktes1,j)- &
                               var(i,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1)+ &
                               var(i-1,ktes1,j)+(var(i-1,ktes1,j)- &
                               var(i-1,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1))
      ENDIF
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts, ktf
   DO i = i_start, i_end + 1
      tmpzx = 0.5*( zx(i,k,j)+ zx(i,k+1,j))
      rdzu = 2./(1./rdzw(i,k,j) + 1./rdzw(i-1,k,j))
      H1(i,k,j)=msfuy(i,j)*(                      &
                 rdx*(var(i,k,j)-var(i-1,k,j)) - tmpzx*         &
                     (H1avg(i,k+1,j)-H1avg(i,k,j))*rdzu )

   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts, ktf
   DO i = i_start, i_end
      dVARdX(i,k,j)= 0.5*(H1(i+1,k,j) + H1(i,k,j))
   ENDDO
   ENDDO
   ENDDO

   END SUBROUTINE calc_dVARdX 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE calc_dVARdY(dVARdY, mu, config_flags, var, iflux,  &
                       msftx, msfty, msfux, msfuy,               &
                       msfvx, msfvy, rdx, rdy,                   &
                       fnm, fnp, cf1, cf2, cf3,                  &
                       zx, zy, rdz, rdzw, dnw, dn,               &
                       ids, ide, jds, jde, kds, kde,             &
                       ims, ime, jms, jme, kms, kme,             &
                       its, ite, jts, jte, kts, kte)

  IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   REAL , INTENT(IN   )           ::        cf1, cf2, cf3

   LOGICAL , INTENT(IN   )  :: iflux  ! .true. : var is a flux
                                      ! .false.: var is no flux

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::     dn
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    dnw

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfux
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfuy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msftx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfty

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   mu
   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: dVARdY
   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::         &
                                                                     rdz, &
                                                                     rdzw

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::    var, &
                                                                     zx, &
                                                                     zy

   REAL ,                                        INTENT(IN   ) ::    rdx, &
                                                                     rdy
! Local data
   INTEGER :: i, j, k, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   REAL , DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1)    ::     H1avg, &
                                                                  H2avg, &
                                                                  H1, &
                                                                  H2
   REAL    :: mrdx, mrdy
   REAL    :: tmpzx, tmpzy, tmpzeta_z, rdzu, rdzv
   INTEGER :: ktes1,ktes2

! End declarations.

!======================================================================

   ktf=MIN(kte,kde-1)

   ktes1=kte-1
   ktes2=kte-2

   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
   IF ( config_flags%periodic_x ) i_start = its
   IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

   ! H2 = partial var over partial y

   DO j = j_start, j_end + 1
   DO k = kts+1,   ktf
   DO i = i_start, i_end
      H2avg(i,k,j)=0.5*(fnm(k)*(var(i,k  ,j-1)+var(i,k  ,j))+  &
                        fnp(k)*(var(i,k-1,j-1)+var(i,k-1,j)))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end + 1
   DO i = i_start, i_end
      IF (iflux .eqv. .true.) THEN
         H2avg(i,kts  ,j)=0.
         H2avg(i,ktf+1,j)=0.
      ELSEIF(iflux .eqv. .false.) THEN
               H2avg(i,kts  ,j)=0.5*(cf1*var(i,1,j  )+cf2*var(i  ,2,j)+ &
                               cf3*var(i,3,j  )+cf1*var(i,1,j-1)+  &
                               cf2*var(i,2,j-1)+cf3*var(i,3,j-1))
         H2avg(i,ktf+1,j)=0.5*(var(i,ktes1,j)+(var(i,ktes1,j)- &
                               var(i,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1)+ &
                               var(i,ktes1,j-1)+(var(i,ktes1,j-1)- &
                               var(i,ktes2,j-1))*0.5*dnw(ktes1)/dn(ktes1))
      ENDIF
   ENDDO
   ENDDO

   DO j = j_start, j_end + 1
   DO k = kts, ktf
   DO i = i_start, i_end
      tmpzy = 0.5*( zy(i,k,j)+ zy(i,k+1,j))
      rdzv = 2./(1./rdzw(i,k,j) + 1./rdzw(i,k,j-1))
      H2(i,k,j)=msfvy(i,j)*(                       &
                 rdy*(var(i,k,j)-var(i,k,j-1)) - tmpzy*          &
                     (H2avg(i ,k+1,j)-H2avg(i,k,j))*rdzv)

   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts, ktf
   DO i = i_start, i_end
      dVARdY(i,k,j)= 0.5*(H2(i,k,j+1) + H2(i,k,j))
   ENDDO
   ENDDO
   ENDDO


   END SUBROUTINE calc_dVARdY 

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE calc_dVARdZ(dVARdZ, config_flags, var, ZP, staggered,    &
                          dn, dnw, rdz, rdzw, fnm, fnp,            &
                          cf1, cf2, cf3,                           &
                          ids, ide, jds, jde, kds, kde,            &
                          ims, ime, jms, jme, kms, kme,            &
                          its, ite, jts, jte, kts, kte)


   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags
   INTEGER ,         INTENT(IN   ) ::       ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte
   REAL , INTENT(IN   )           ::        cf1, cf2, cf3
   LOGICAL, INTENT(IN) :: staggered  ! true:  staggered    => output nonstaggered
                                     ! false: nonstaggered => output staggered
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp
   REAL , DIMENSION( kms:kme ) ,            INTENT(IN   )      ::  dn
   REAL , DIMENSION( kms:kme ) ,            INTENT(IN   )      :: dnw

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: dVARdZ

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme ) ,                       &
                                            INTENT(IN   )      ::   var, &
                                                                    rdz, &
                                                                    rdzw, &
                                                                    ZP
! LOCAL VARS
   INTEGER :: i, j, k, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   REAL , DIMENSION( its:ite, kts:kte, jts:jte)            ::        H3
! End declarations.
!-----------------------------------------------------------------------
   ktf=MIN(kte,kde-1)
   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

! H3
   IF (staggered .eqv. .true.) THEN
      DO j = j_start, j_end
      DO k = kts,ktf
      DO i = i_start, i_end
         H3(i,k,j)=(var(i,k+1,j)-var(i,k,j))*rdzw(i,k,j)
      ENDDO
      ENDDO
      ENDDO
   ENDIF

   IF (staggered .eqv. .false.) THEN
      DO j = j_start, j_end
      DO k = kts+1,ktf
      DO i = i_start, i_end
         H3(i,k,j)=(var(i,k,j)-var(i,k-1,j))*rdz(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j = j_start, j_end
      DO i = i_start, i_end
!         H3(i,kts,j)=cf1*H3(i,2,j)+cf2*H3(i,3,j)+cf3*H3(i,4,j)
!         H3(i,kts,j)=0.
         H3(i,kts,j)=H3(i,kts+1,j) ! Leukauf

!         z0 = ZPW(i,1,j)
!         z1 = ZPM(i,1,j)
!         z2 = ZPM(i,2,j)
!         w1 = (z0 - z2)/(z1 - z2)
!         w2 = 1. - w1
!         var0 = w1*var(i,1,j)+w2*var(i,2,j)
!         H3(i,kts,j) = (var(i,1,j)-var0*rdzw(i,k,j)
 
!         H3(i,kts,j)=(var(i,kts,j)-varsfc(i,j))*2.*rdzw(i,kts,j)
!         H3(i,kts,j)=cf1*var(i  ,1,j)+cf2*var(i  ,2,j)+cf3*var(i  ,3,j)
!         H3(i,kts,j)=(var(i,kts,j)-H3(i,kts,j))*2.*rdzw(i,kts,j)
!          H3(i,kts,j) = interp_0( H3(i,kts+1:ktf,j), ZP(i,kts+1:ktf,j),
!          ZP(i,kts,j), kme-kms+1 ) 
         H3(i,ktf+1,j)=0.
      ENDDO
      ENDDO
   ENDIF

   DO j = j_start, j_end
   DO k = kts,ktf
   DO i = i_start, i_end
      dVARdZ(i,k,j)=H3(i,k,j)
   ENDDO
   ENDDO
   ENDDO

   END SUBROUTINE calc_dVARdZ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE horizontal_diffusion_sX(tendency, config_flags, var, flx,      &
                                   msftx, msfty, msfux, msfuy,            &
                                   msfvx, msfvy, xkhh, rdx, rdy,          &
                                   fnm, fnp, cf1, cf2, cf3,               &
                                   zx, zy, rdz, rdzw, dnw, dn,            &
                                   ids, ide, jds, jde, kds, kde,          &
                                   ims, ime, jms, jme, kms, kme,          &
                                   its, ite, jts, jte, kts, kte)

!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

   REAL , INTENT(IN   )           ::        cf1, cf2, cf3

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::     dn
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    dnw

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfux
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfuy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msftx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfty


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: tendency, flx

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::         &
                                                                    xkhh, &
                                                                     rdz, &
                                                                     rdzw

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::    var, &
                                                                      zx, &
                                                                      zy

   REAL ,                                        INTENT(IN   ) ::    rdx, &
                                                                     rdy

! Local data

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1)    ::     H1avg, &
                                                                     H1, &
                                                                 xkxavg

   REAL    :: mrdx, mrdy, rcoup
   REAL    :: tmpzx, tmpzy, tmpzeta_z, rdzu, rdzv
   INTEGER :: ktes1,ktes2

! End declarations.
!-----------------------------------------------------------------------

   ktf=MIN(kte,kde-1)
 
!-----------------------------------------------------------------------
! scalars:   t (.), u(|), v(+), w(-)
!       
!       t  u  t  u                               t  v  t  v 
!
! w     -     3     -      k+1             w     -     3     -      k+1 
!
! t     .  1  O  1  .      k               t     .  2  O  2  .      k      
!
! w     -     3     -      k               w     -     3     -      k   
!
! t     .  |  .  |  .      k-1             t     .  +  .  +  .      k-1 
!
! w     -     -     -      k-1             w     -     -     -      k-1 
!
! t    i-1 i  i i+1                             j-1 j  j j+1         
!

   ktes1=kte-1
   ktes2=kte-2

   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

! H1 = partial var over partial x

   DO j = j_start, j_end
   DO k = kts, ktf
   DO i = i_start, i_end + 1
      xkxavg(i,k,j)=0.5*(xkhh(i-1,k,j)+xkhh(i,k,j))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts+1, ktf
   DO i = i_start, i_end + 1
      H1avg(i,k,j)=0.5*(fnm(k)*(var(i-1,k  ,j)+var(i,k  ,j))+  &
                        fnp(k)*(var(i-1,k-1,j)+var(i,k-1,j)))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO i = i_start, i_end + 1
      H1avg(i,kts  ,j)=0.5*(cf1*var(i  ,1,j)+cf2*var(i  ,2,j)+ &
                            cf3*var(i  ,3,j)+cf1*var(i-1,1,j)+  &
                            cf2*var(i-1,2,j)+cf3*var(i-1,3,j))
      H1avg(i,ktf+1,j)=0.5*(var(i,ktes1,j)+(var(i,ktes1,j)- &
                            var(i,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1)+ &
                            var(i-1,ktes1,j)+(var(i-1,ktes1,j)- &
                            var(i-1,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1))
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts, ktf
   DO i = i_start, i_end + 1
! new
      tmpzx = 0.5*( zx(i,k,j)+ zx(i,k+1,j))
      rdzu = 2./(1./rdzw(i,k,j) + 1./rdzw(i-1,k,j))
      H1(i,k,j)=-msfuy(i,j)*xkxavg(i,k,j)*(                      &
                 rdx*(var(i,k,j)-var(i-1,k,j)) - tmpzx*         &
                     (H1avg(i,k+1,j)-H1avg(i,k,j))*rdzu )

   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts+1, ktf
   DO i = i_start, i_end
      H1avg(i,k,j)=0.5*(fnm(k)*(H1(i+1,k  ,j)+H1(i,k  ,j))+  &
                        fnp(k)*(H1(i+1,k-1,j)+H1(i,k-1,j)))

      tmpzx = 0.5*( zx(i,k,j)+ zx(i+1,k,j  ))

      H1avg(i,k,j)=H1avg(i,k,j)*tmpzx

   ENDDO
   ENDDO
   ENDDO
 
   DO j = j_start, j_end
   DO i = i_start, i_end
      H1avg(i,kts  ,j)=0.
      H1avg(i,ktf+1,j)=0.
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts,ktf
   DO i = i_start, i_end

      mrdx=msftx(i,j)*rdx
      mrdy=msfty(i,j)*rdy
      
      flx(i,k,j) = H1(i,k,j)

      tendency(i,k,j)=                                           &
          -(mrdx*0.5*(H1(i+1,k,j)-H1(i  ,k,j))-                  &
           msfty(i,j)*(H1avg(i,k+1,j)-H1avg(i,k,j))*rdzw(i,k,j))

   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE horizontal_diffusion_sX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE horizontal_diffusion_sY(tendency, config_flags, var, flx,      &
                                   msftx, msfty, msfux, msfuy,            &
                                   msfvx, msfvy, xkhh, rdx, rdy,          &
                                   fnm, fnp, cf1, cf2, cf3,               &
                                   zx, zy, rdz, rdzw, dnw, dn,            &
                                   ids, ide, jds, jde, kds, kde,          &
                                   ims, ime, jms, jme, kms, kme,          &
                                   its, ite, jts, jte, kts, kte)

!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                           ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

   REAL , INTENT(IN   )           ::        cf1, cf2, cf3

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::     dn
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    dnw

   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfux
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfuy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfvy
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msftx
   REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN   ) ::   msfty

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: tendency, flx

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::         &
                                                                    xkhh, &
                                                                     rdz, &
                                                                     rdzw

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::    var, &
                                                                      zx, &
                                                                      zy

   REAL ,                                        INTENT(IN   ) ::    rdx, &
                                                                     rdy

! Local data

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1)    ::     H2avg, &
                                                                     H2, &
                                                                 xkxavg

   REAL    :: mrdx, mrdy, rcoup
   REAL    :: tmpzx, tmpzy, tmpzeta_z, rdzu, rdzv
   INTEGER :: ktes1,ktes2

! End declarations.
!-----------------------------------------------------------------------

   ktf=MIN(kte,kde-1)
 
!-----------------------------------------------------------------------
! scalars:   t (.), u(|), v(+), w(-)
!       
!       t  u  t  u                               t  v  t  v 
!
! w     -     3     -      k+1             w     -     3     -      k+1 
!
! t     .  1  O  1  .      k               t     .  2  O  2  .      k      
!
! w     -     3     -      k               w     -     3     -      k   
!
! t     .  |  .  |  .      k-1             t     .  +  .  +  .      k-1 
!
! w     -     -     -      k-1             w     -     -     -      k-1 
!
! t    i-1 i  i i+1                             j-1 j  j j+1         
!

   ktes1=kte-1
   ktes2=kte-2

   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

! H2 = partial var over partial y

   DO j = j_start, j_end + 1
   DO k = kts, ktf
   DO i = i_start, i_end
      xkxavg(i,k,j)=0.5*(xkhh(i,k,j-1)+xkhh(i,k,j))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end + 1
   DO k = kts+1,   ktf
   DO i = i_start, i_end
      H2avg(i,k,j)=0.5*(fnm(k)*(var(i,k  ,j-1)+var(i,k  ,j))+  &
                        fnp(k)*(var(i,k-1,j-1)+var(i,k-1,j)))
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end + 1
   DO i = i_start, i_end
      H2avg(i,kts  ,j)=0.5*(cf1*var(i,1,j  )+cf2*var(i  ,2,j)+ &
                            cf3*var(i,3,j  )+cf1*var(i,1,j-1)+  &
                            cf2*var(i,2,j-1)+cf3*var(i,3,j-1))
      H2avg(i,ktf+1,j)=0.5*(var(i,ktes1,j)+(var(i,ktes1,j)- &
                            var(i,ktes2,j))*0.5*dnw(ktes1)/dn(ktes1)+ &
                            var(i,ktes1,j-1)+(var(i,ktes1,j-1)- &
                            var(i,ktes2,j-1))*0.5*dnw(ktes1)/dn(ktes1))
   ENDDO
   ENDDO

   DO j = j_start, j_end + 1
   DO k = kts, ktf
   DO i = i_start, i_end
      tmpzy = 0.5*( zy(i,k,j)+ zy(i,k+1,j))
      rdzv = 2./(1./rdzw(i,k,j) + 1./rdzw(i,k,j-1))
      H2(i,k,j)=-msfvy(i,j)*xkxavg(i,k,j)*(                       &
                 rdy*(var(i,k,j)-var(i,k,j-1)) - tmpzy*          &
                     (H2avg(i ,k+1,j)-H2avg(i,k,j))*rdzv)

   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts+1, ktf
   DO i = i_start, i_end
      H2avg(i,k,j)=0.5*(fnm(k)*(H2(i,k  ,j+1)+H2(i,k  ,j))+  &
                        fnp(k)*(H2(i,k-1,j+1)+H2(i,k-1,j)))

      tmpzy = 0.5*( zy(i,k,j)+ zy(i  ,k,j+1))

      H2avg(i,k,j)=H2avg(i,k,j)*tmpzy

   ENDDO
   ENDDO
   ENDDO
 
   DO j = j_start, j_end
   DO i = i_start, i_end
      H2avg(i,kts  ,j)=0.
      H2avg(i,ktf+1,j)=0.
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts,ktf
   DO i = i_start, i_end

      mrdx=msftx(i,j)*rdx
      mrdy=msfty(i,j)*rdy

      flx(i,k,j) = H2(i,k,j)

      tendency(i,k,j)=                                      &
          -(mrdy*0.5*(H2(i,k,j+1)-H2(i,k,j  ))-             &
           msfty(i,j)*(H2avg(i,k+1,j)-H2avg(i,k,j))*rdzw(i,k,j))
   
   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE horizontal_diffusion_sY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE vertical_diffusion_sZ( tendency, config_flags, var, flx, xkhv, sfcflux, &
                                 rdz, rdzw, fnm, fnp,                         &
                                 ids, ide, jds, jde, kds, kde,                &
                                 ims, ime, jms, jme, kms, kme,                &
                                 its, ite, jts, jte, kts, kte)
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,         INTENT(IN   ) ::       ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnm
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::    fnp

   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::tendency, flx

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(IN) ::   xkhv

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::   sfcflux

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme ) ,                       &
                                            INTENT(IN   )      ::    var, &
                                                                     rdz, &
                                                                     rdzw

! LOCAL VARS

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION( its:ite, kts:kte, jts:jte)            ::        H3, &
                                                                 xkxavg, &
                                                                  rravg

! End declarations.
!-----------------------------------------------------------------------


   ktf=MIN(kte,kde-1)
  
   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-2,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

! H3

   xkxavg = 0.

   DO j = j_start, j_end
   DO k = kts+1,ktf
   DO i = i_start, i_end
      xkxavg(i,k,j)=fnm(k)*xkhv(i,k,j)+fnp(k)*xkhv(i,k-1,j)
      H3(i,k,j)=-xkxavg(i,k,j)*(var(i,k,j)-var(i,k-1,j))*rdz(i,k,j)
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO i = i_start, i_end
      H3(i,kts,j)=0.
      H3(i,ktf+1,j)=0.
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO k = kts,ktf
   DO i = i_start, i_end
      flx(i,k,j) = H3(i,k,j)

      tendency(i,k,j)= -(H3(i,k+1,j)-H3(i,k,j))*rdzw(i,k,j)
   ENDDO
   ENDDO
   ENDDO

   DO j = j_start, j_end
   DO i = i_start, i_end
      flx(i,kts,j) = sfcflux(i,j)

      tendency(i,kts,j)= tendency(i,kts,j) + sfcflux(i,j)*rdzw(i,kts,j)
   ENDDO
   ENDDO

END SUBROUTINE vertical_diffusion_sZ

end MODULE module_avgflx_em

