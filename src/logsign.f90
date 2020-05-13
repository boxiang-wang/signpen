! --------------------------------------------------------------------------
! log_sign.f90: the perceptron penalty for sign consistency of logistic regression. 
! --------------------------------------------------------------------------
!
! USAGE:
! 
! call log_sign (lam2, lam3, gam, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
! & pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, &
! & ibeta, nbeta, alam, npass, jerr) 
! 
! INPUT ARGUMENTS:
! 
!    lam2 = regularization parameter for the quadratic penalty 
!    lam3 = regularization parameter for the perceptron penalty
!    gam(nvars) = prior information on the signs of beta
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor {-1, 1} 
!            for classification.
!    jd(jd(1)+1) = predictor variable deletion flag
!                  jd(1) = 0  => use all variables
!                  jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                pf2(j) = 0 => jth variable unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original
!                variables locations and scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nalam) = intercept values for each solution
!    beta(pmax, nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)

SUBROUTINE log_sign (lam2, lam3, gam, nobs, nvars, x, y, jd, pf, pf2, & 
   & dfmax, pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, &
   & b0, beta, ibeta, nbeta, alam, npass, jerr) 

   IMPLICIT NONE
   !------- arg types ----------------------------------------------
   INTEGER :: nobs, nvars, nlam, nalam, maxit, dfmax, pmax
   INTEGER :: isd, npass, jerr, jd (*), ibeta (pmax), nbeta(nlam)
   DOUBLE PRECISION :: lam2, lam3, flmin, eps
   DOUBLE PRECISION :: x (nobs, nvars), y (nobs)
   DOUBLE PRECISION :: gam (nvars)
   DOUBLE PRECISION :: pf (nvars), pf2 (nvars)
   DOUBLE PRECISION :: ulam (nlam), alam (nlam)
   DOUBLE PRECISION :: beta (pmax, nlam), b0 (nlam)


   !------- local declarations -------------------------------------
   INTEGER :: j, l, nk, ju (nvars)
   DOUBLE PRECISION :: xmean (nvars), xnorm (nvars), maj (nvars)
   LOGICAL :: strong

   !------- preliminary step ---------------------------------------
   CALL chkvars (nobs, nvars, x, ju)
   IF (jd(1) > 0) ju(jd(2:(jd(1)+1))) = 0
   IF (Maxval (ju) <= 0) THEN
      jerr = 7777
      RETURN
   ENDIF
   IF (Maxval (pf) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF
   IF (Maxval (pf2) <= 0.0D0) THEN
      jerr = 10000
      RETURN
   ENDIF
   pf = Max (0.0D0, pf)
   pf2 = Max (0.0D0, pf2)
   strong = .TRUE.

   !------- first standardize the data -----------------------------
   CALL Standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
   CALL log_sign_path (lam2, lam3, gam, maj, nobs, nvars, x, y, ju, &
      & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, strong, &
      & nalam, b0, beta, ibeta, nbeta, alam, npass, jerr)
   IF (jerr > 0) RETURN  ! check error after calling function

   !------- organize beta afterward --------------------------------
   DO l = 1, nalam
      nk = nbeta(l)
      IF (isd == 1) THEN
         DO j = 1, nk
            beta(j, l) = beta(j, l) / xnorm(ibeta(j))
         ENDDO
      ENDIF
      b0(l) = b0(l) - Dot_product (beta(1:nk, l), &
       & xmean(ibeta(1:nk)))
   ENDDO
   RETURN
END SUBROUTINE log_sign

SUBROUTINE log_sign_path (lam2, lam3, gam, maj, nobs, nvars, x, y, &
   & ju, pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, &
   & strong, nalam, b0, beta, m, nbeta, alam, npass, jerr)

   IMPLICIT NONE
   !------- arg types -----------------------------------------------
   DOUBLE PRECISION, PARAMETER :: BIG = 9.9E30
   DOUBLE PRECISION, PARAMETER :: MFL = 1.0E-6
   INTEGER, PARAMETER :: MNLAM = 6
   INTEGER :: nobs, nvars, dfmax, pmax, nlam
   INTEGER :: maxit, nalam, npass, jerr 
   INTEGER :: ju (nvars), m (pmax), nbeta (nlam)
   DOUBLE PRECISION :: lam2, lam3, flmin, eps
   DOUBLE PRECISION :: x (nobs, nvars), y (nobs)
   DOUBLE PRECISION :: pf (nvars), pf2 (nvars), gam (nvars), ulam (nlam)
   DOUBLE PRECISION :: beta (pmax, nlam), b0 (nlam)
   DOUBLE PRECISION :: alam (nlam), maj (nvars)
   LOGICAL :: strong

   !------- local declarations -------------------------------------
   DOUBLE PRECISION :: d, dif, oldb, u, v, al, alf
   DOUBLE PRECISION :: dl (nobs), r (nobs)
   DOUBLE PRECISION :: b (0:nvars), oldbeta (0:nvars)
   INTEGER :: i, k, j, l, ctr, ni, me, mnl, mm (nvars)

   !------- local declarations for the strong rule -----------------
   INTEGER :: jx, jxx (nvars)
   DOUBLE PRECISION :: tlam, ga (nvars), vl (nvars), al0

   !------- some initial setup ------------------------------------- 
   al = 0.0D0
   r = 0.0D0
   b = 0.0D0
   oldbeta = 0.0D0
   m = 0
   mm = 0
   npass = 0
   ni = npass
   maj = maj * 0.25
   alf = 0.01D0

   ! strong rule
   !   jxx = 0 using the strong rule
   !   jxx = 1 not using the strong rule
   jxx = 0
   IF (strong .EQV. .TRUE.) THEN
      jxx = 0
   ELSE
      jxx = 1
   ENDIF

   !---------- lambda loop -----------------------------------------
   mnl = Min (MNLAM, nlam)
   IF (flmin < 1.0D0) THEN
      flmin = Max (MFL, flmin)
      alf = flmin ** (1.0D0 / (nlam-1.0D0))
   ENDIF
   vl = 0.0D0
   vl = dot_product(y/(1.0D0+Exp(r)), x(:, j))
   !CALL logdrv(nobs, nvars, x, y, r, vl)
   ga = Abs (vl)

   loop_lambda: DO l = 1, nlam
      al0 = al
      IF (flmin >= 1.0D0) THEN
         al = ulam(l)
      ELSE
         IF (l > 2) THEN
            al = al * alf
         ELSE IF (l == 1) THEN
            al = BIG
         ELSE IF (l == 2) THEN
            al0 = 0.0D0
            DO j = 1, nvars
               IF (ju(j) /= 0) THEN
                  IF (pf(j) > 0.0D0) THEN
                     al0 = Max (al0, Abs(dot_product(y/(1.0D0+Exp(r)), &
                      & x(:, j)))/pf(j))
                  ENDIF
               ENDIF
            ENDDO
            al = al0 * alf
         ENDIF
      ENDIF
      ctr = 0

      !---------- check strong rule (in lamdba loop) ------------------
      tlam = 2.0 * al - al0
      loop_strong_rule: DO j = 1, nvars
         IF (jxx(j) == 1) CYCLE
         IF (ga(j) > pf(j) * tlam) jxx(j) = 1
      ENDDO loop_strong_rule

      !---------- outer loop ------------------------------------------
      loop_outer: DO
         oldbeta(0) = b(0)
         IF (ni > 0) oldbeta (m(1:ni)) = b(m(1:ni))

         !---------- middle loop -----------------------------------------
         loop_middle: DO
            npass = npass + 1
            dif = 0.0D0

            !---------- 1. update beta_j first ( in middle loop) ------------
            loop_middle_update_betaj: DO k = 1, nvars
               IF (jxx(k) == 0) CYCLE
               IF (ju(k) /= 0) THEN
                  oldb = b(k)
                  u = 0.0D0
                  ! Note that most literature define the margin as u, 
                  !    but the margin is actually r in the code.
                  ! u: 
                  ! \sum_{i=1}^n V'(r_i)y_i x_{ij}
                  DO i = 1, nobs   
                     dl(i) = -1.0D0 / (1.0D0 + Exp(r(i)))
                     u = u + dl(i) * y(i) * x(i, k)
                  ENDDO 
                  ! u:
                  ! M \tilde{\beta_j} - 1/n * (\sum_{i=1}^n V'(r_i)y_i x_{ij} )
                  u = maj(k) * b(k) - u / nobs
                  v = al * pf(k)
                  IF (gam(k) * u < 0) v = v + lam3 * Abs(gam(k))
                  v = Abs (u) - v
                  ! The update of the coefficients in Algorithm 1:
                  IF (v > 0.0D0) THEN
                     b(k) = Sign (v, u) / (maj(k) + pf2(k) * lam2)
                  ELSE
                     b(k) = 0.0D0
                  ENDIF
                  d = b(k) - oldb
                  IF (Abs (d) > 0.0D0) THEN
                     dif = Max (dif, d * d)  !!!
                     r = r + y * x(:, k) * d
                     IF (mm(k) == 0) THEN
                        ni = ni + 1
                        IF (ni > pmax) EXIT
                        mm(k) = ni
                        m(ni) = k !indicate which one is non-zero
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO loop_middle_update_betaj

            !---------- 2. update intercept ( in middle loop) ---------------
            IF (ni > pmax) EXIT
            d = 0.0D0
            DO i = 1, nobs
               dl(i) = -1.0D0 / (1.0D0 + Exp(r(i)))
               d = d + dl(i) * y(i)
            ENDDO
            d = - 4.0D0 * d / nobs !!!!!!!!!
            IF (d /= 0.0D0) THEN
               ! New intercept:
               ! newbeta0 = oldbeta0 - 1/M * \sum_{i=1}^n V'(r_i) y_i / n                  
               b(0) = b(0) +  d
               r = r + y * d
               dif = Max (dif, d * d) !!!
            ENDIF
            IF (dif < eps) EXIT
            IF (npass > maxit) THEN
               jerr = -1
               RETURN
            ENDIF

            !---------- inner loop ------------------------------------------
            loop_inner: DO
               npass = npass + 1
               dif = 0.0D0

               !---------- 1. update beta_j first ( in inner loop) -------------
               DO j = 1, ni
                  k = m(j)
                  oldb = b(k)
                  u = 0.0D0
                  DO i = 1, nobs   !!!!!!!!!!!!!!!!!!!!!!
                     dl(i) = -1.0D0 / (1.0D0 + Exp(r(i)))
                     u = u + dl(i) * y(i) * x(i, k)
                  ENDDO !!!!!!!!!!!!!!!!!!!!!!
                  u = maj(k) * b(k) - u / nobs
                  v = al * pf(k)
                  IF (gam(k) * u < 0) v = v + lam3 * Abs(gam(k))
                  v = Abs (u) - v
                  IF (v > 0.0D0) THEN
                     b(k) = Sign (v, u) / (maj(k) + pf2(k) * lam2)
                  ELSE
                     b(k) = 0.0D0
                  ENDIF
                  d = b(k) - oldb
                  IF (Abs(d) > 0.0D0) THEN
                     dif = Max (dif, d * d) !!!
                     r = r + y * x(:, k) * d
                  ENDIF
               ENDDO
               !---------- 2. update intercept ( in inner loop) ----------------
               d = 0.0D0
               DO i = 1, nobs   !!!!!!!!!!!!!!!!!!!!!!
                  dl(i) = -1.0D0 / (1.0D0 + Exp(r(i)))
                  d = d + dl(i) * y(i)
               ENDDO !!!!!!!!!!!!!!!!!!!!!!
               d = - 4.0D0 * d / nobs !!!!!!!!!
               IF (d /= 0.0D0) THEN
                  b(0) = b(0) + d
                  r = r + y * d
                  dif = Max (dif, d * d) !!!
               ENDIF
               IF (dif < eps) EXIT
               IF (npass > maxit) THEN
                  jerr = -1
                  RETURN
               ENDIF
            ENDDO loop_inner
         ENDDO loop_middle
         IF (ni > pmax) EXIT

         !---------- final check -----------------------------------------
         jx = 0
         IF ((b(0) - oldbeta(0)) ** 2 >= eps) jx = 1
         DO j = 1, ni
            IF ((b(m(j)) - oldbeta(m(j))) ** 2 >= eps) THEN
               jx = 1
               EXIT
            ENDIF
         ENDDO
         IF (jx /= 0) CYCLE
         !---------- check the KKT conditions of the discarded variables--    
         CALL logdrv (nobs, nvars, x, y, r, vl)
         ga = Abs(vl)
         DO j = 1, nvars
            IF (jxx(j) == 1) CYCLE
            If (ga(j) > al * pf(j)) THEN
               jxx(j) = 1
               jx = 1
            ENDIF
         ENDDO
         ! jx:
         ! jx checks if some variable violates the KKT conditions. 
         IF (jx == 1) CYCLE
         EXIT
         ctr = ctr + 1
         IF (ctr > maxit) THEN
            jerr = - l
            RETURN
         ENDIF
      ENDDO loop_outer

      !---------- final update variable save results ------------------
      IF (ni > pmax) THEN
         jerr = - 10000 - l
         EXIT
      ENDIF
      IF (ni > 0) beta(1:ni, l) = b(m(1:ni))
      nbeta(l) = ni
      b0(l) = b(0)
      alam(l) = al
      nalam = l
      IF (l < mnl) CYCLE
      IF (flmin >= 1.0D0) CYCLE
      me = Count (beta(1:ni, l) /= 0.0D0)
      IF (me > dfmax) EXIT
   ENDDO loop_lambda
   RETURN
END SUBROUTINE log_sign_path

! SUBROUTINE logdrv (nobs, nvars, x, y, r, vl)
!    IMPLICIT NONE
!    INTEGER :: nobs, nvars, i
!    DOUBLE PRECISION :: x (nobs, nvars), y (nobs)
!    DOUBLE PRECISION :: r (nobs), vl (nvars), dl (nobs), dly (nobs)
!    vl = 0.0D0
!    dl = 0.0D0
!    DO i = 1, nobs
!       dl(i) = -1.0D0 / (1.0D0 + Exp(r(i))
!    ENDDO
!    dly = dl * y
!    DO i = 1, nvars
!       vl(i) = Dot_product(dly, x(1:nobs, i))
!    ENDDO
! END SUBROUTINE logdrv

