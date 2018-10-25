subroutine ausmplusup()
  !
  ! Implementation of AUSM+Up
  !

  use m_init

  implicit none

  integer :: j
  real( kind=dp ) :: rhol, rhor, ul, ur, htl, htr
  real( kind=dp ) :: astarsqrl, astarsqrr, astarl, astarr, atildel, atilder, ahalf
  real( kind=dp ) :: Mj, Mjp1, Mplus, Mminus, Mhalf, mjplushalfplus, mjplushalfminus
  real( kind=dp ) :: pl, pr, pj, pjp1, Pplus, Pminus, Phalf
  real( kind=dp ) :: Mzero, Mzerosqr, fa, ahatl, ahatr, rhohalf
  real( kind=dp ) :: alpha, sigma, ku, kp, mdothalf, mbarsqr, mhatsqr
  real( kind=dp ) :: v1, p1, a1, Minf, Minfsqr
  
  do j = 0, ncells+1

     ! initializations ---------------------------------------------------------

     sigma = 1.
     kp    = 0.25
     ku    = 0.75

     rhol = u1 (j)
     rhor = u1 (j+1)
     
     ul  = u2 (j) / rhol
     ur  = u2 (j+1) / rhor

     pl = (gamma-1.) * ( u3 (j) - 0.5 * u2 (j) * u2 (j) / u1 (j) )
     pr = (gamma-1.) * ( u3 (j+1) - 0.5 * u2 (j+1) * u2 (j+1) / u1 (j+1) )   

     htl = (u3 (j) + pl) / rhol
     htr = (u3 (j+1) + pr) / rhor

     v1 = u2 (1) / u1 (1)
     p1 = (gamma-1.) * ( u3 (1) - 0.5 * u2 (1) * u2 (1) / u1 (1) )
     a1 = sqrt( ( gamma * p1 ) / u1 (1) )
     Minf = v1 / a1
     Minfsqr = Minf * Minf

     ! compute interface speed of sound ----------------------------------------
     
     astarsqrl = 2. * (gamma - 1.) * htl * (gamma + 1.)     ! (29) [L06]
     astarsqrr = 2. * (gamma - 1.) * htr * (gamma + 1.)

     astarl = sqrt( astarsqrl ) ! (28) [L06]
     astarr = sqrt( astarsqrr )

     ahatl = astarsqrl / ( max( astarl, abs( ul ) ) )       ! (28) [L06]
     ahatr = astarsqrr / ( max( astarr, abs( ur ) ) )

     ahalf = min( ahatl, ahatr )
     
     ! compute interface Mach number -------------------------------------------

     Mj   = ul / ahalf
     Mjp1 = ur / ahalf

     Mbarsqr = ( ul*ul + ur*ur ) / ( 2. * ahalf * ahalf)    ! (70) [L06]
     Mzerosqr = min( 1., max( Mhatsqr, Minfsqr ) )          ! (71) [L06]
     Mzero = sqrt( Mzerosqr )
     fa = Mzero * ( 2. - Mzero )                            ! (72) [L06]

     if (abs( Mj) >= 1.) then
        Mplus = 0.5* ( Mj + abs( Mj) )
     else
        Mplus = 0.25 * (Mj +1.)*(Mj +1.) &
             + 0.125 * (Mj*Mj-1.)* (Mj*Mj-1.) ! 1/4 (19) [L06]
     end if

     if (abs( Mjp1 ) >= 1.) then
        Mminus = 0.5* ( Mjp1 + abs( Mjp1 ) )
     else
        Mminus = -0.25 * (Mjp1 - 1.) * (Mjp1 - 1.) - 0.125 &
             * (Mjp1*Mjp1 - 1.)* (Mjp1*Mjp1 - 1.)
     end if

     rhohalf = (rhol + rhor ) / 2.

     Mhalf = Mplus + Mminus - (kp / fa) * max( 1.-sigma*Mhatsqr, 0.) &
          * (pr - pl)/(rhohalf * ahalf * ahalf) ! (73) [L06]

     ! compute mass flux -------------------------------------------------------

     if (Mhalf > 0.) then
        mdothalf = ahalf * Mhalf * rhol
     else
        mdothalf = ahalf * Mhalf * rhor
     end if
     
     ! compute pressure flux ---------------------------------------------------

     pj   = pl
     pjp1 = pr
     alpha = 5. * fa * fa - 4.

     if (abs( Mj ) >= 1.) then
        Pplus = 0.5* ( 1. + sign( 1.0_dp, Mj) )
     else
        Pplus = 0.25 * (Mj +1.) * (Mj +1.) * ( 2. - Mj) &
             + 0.1875 * alpha * Mj * (Mj*Mj - 1.) * (Mj*Mj - 1.)
     end if

     if (abs( Mjp1 ) >= 1.) then
        Pminus = 0.5* ( 1. - sign( 1.0_dp, Mjp1) )
     else
        Pminus = 0.25 * (Mjp1 - 1.)*(Mjp1 - 1.) * ( 2. + Mjp1) &
             - 0.1875 * alpha * Mjp1 * (Mjp1*Mjp1 - 1.) * (Mjp1*Mjp1 - 1.)
     end if     
     
     Phalf = Pplus * pj + Pminus * pjp1 - ku * Pplus * Pminus &
          * (rhol + rhor) * (fa*ahalf) * (ur-ul)

     ! compute numerical fluxes ------------------------------------------------

     if (mdothalf > 0.) then
        f1 (j) = mdothalf * rhol
        f2 (j) = mdothalf * rhol * ul + Phalf
        f3 (j) = mdothalf * rhol * htl
     else
        f1 (j) = mdothalf * rhor
        f2 (j) = mdothalf * rhor * ur + Phalf
        f3 (j) = mdothalf * rhor * htr        
     end if

  end do

end subroutine ausmplusup
