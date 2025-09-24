SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
               PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, &
               TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, &
               ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, &
               PNEWDT, JPROPS, NJPROP, PERIOD)
      INCLUDE 'ABA_PARAM.INC'

!     --- Inputs / Outputs ---
      INTEGER NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JELEM, NPREDF, MLVARX, NDLOAD, NJPROP
      DOUBLE PRECISION RHS(NDOFEL,NRHS), AMATRX(NDOFEL,NDOFEL)
      DOUBLE PRECISION SVARS(NSVARS), ENERGY(8)
      DOUBLE PRECISION PROPS(NPROPS), COORDS(MCRD,NNODE)
      DOUBLE PRECISION U(NDOFEL), DU(MLVARX,NRHS), V(NDOFEL), A(NDOFEL)
      DOUBLE PRECISION TIME(2), DTIME, PNEWDT, ADLMAG(NDOFEL), DDLMAG(NDOFEL)
      DOUBLE PRECISION PREDEF(NPREDF), MDLOAD(NDLOAD), PARAMS(*)
      INTEGER JTYPE, JDLTYP, PERIOD
      DOUBLE PRECISION JPROPS(NJPROP)

!     --- Local variables ---
      INTEGER i,j,k,l,INTPT
      DOUBLE PRECISION X(NNODE), Y(NNODE)
      DOUBLE PRECISION dN_dxi(NNODE,2), dN_dx(NNODE,2)
      DOUBLE PRECISION JACOB(2,2), JACOB_INV(2,2), CO_FACTOR(2,2), JACOB_DET
      DOUBLE PRECISION BMAT_U(3,2*NNODE)
      DOUBLE PRECISION eps(3,1), STRESS(3), DDSDDE(3,3)
      DOUBLE PRECISION Fint_u(2*NNODE)
      DOUBLE PRECISION gp(9,2), gw(9)
      DOUBLE PRECISION weight
      DOUBLE PRECISION a


!     --- Initialize ---
      AMATRX(:,:) = 0.0D0
      RHS(:,:)   = 0.0D0
      Fint_u(:)  = 0.0D0

!     --- Extract coordinates ---
      DO i = 1, NNODE
          X(i) = COORDS(1,i)
          Y(i) = COORDS(2,i)
      END DO

!     --- 3x3 Gauss points for Q8 element ---
            a = SQRT(3.0D0/5.0D0)
      gp(1,:) = (/ -a, -a /)
      gp(2,:) = (/  0.0D0, -a /)
      gp(3,:) = (/  a, -a /)
      gp(4,:) = (/ -a, 0.0D0 /)
      gp(5,:) = (/  0.0D0, 0.0D0 /)
      gp(6,:) = (/  a, 0.0D0 /)
      gp(7,:) = (/ -a, a /)
      gp(8,:) = (/  0.0D0, a /)
      gp(9,:) = (/  a, a /)

      gw(:) = (/ 5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0, &
                 8.0D0/9.0D0, 16.0D0/9.0D0, 8.0D0/9.0D0, &
                 5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0 /)

!     --- Loop over Gauss points ---
      DO INTPT = 1, 9

!        --- Compute shape derivatives at this Gauss point ---
          CALL SHAPE_FUNCTN(gp(INTPT,1), eta=gp(INTPT,2), dN_dxi)

!        --- Jacobian ---
          JACOB(:,:) = 0.0D0
          DO i=1,NNODE
              JACOB(1,1) = JACOB(1,1) + dN_dxi(i,1)*X(i)
              JACOB(1,2) = JACOB(1,2) + dN_dxi(i,2)*X(i)
              JACOB(2,1) = JACOB(2,1) + dN_dxi(i,1)*Y(i)
              JACOB(2,2) = JACOB(2,2) + dN_dxi(i,2)*Y(i)
          END DO

          JACOB_DET = JACOB(1,1)*JACOB(2,2)-JACOB(1,2)*JACOB(2,1)
          IF (JACOB_DET <= 0.0D0) THEN
              WRITE(*,*) 'Error: Jacobian <=0 at elem', JELEM, 'GP', INTPT
              STOP
          END IF

!         --- Inverse Jacobian ---
          CO_FACTOR(1,1) =  JACOB(2,2)
          CO_FACTOR(2,2) =  JACOB(1,1)
          CO_FACTOR(1,2) = -JACOB(1,2)
          CO_FACTOR(2,1) = -JACOB(2,1)
          DO i=1,2
              DO j=1,2
                  JACOB_INV(i,j) = CO_FACTOR(i,j)/JACOB_DET
              END DO
          END DO

!        --- Compute dN/dx ---
          DO i=1,NNODE
              dN_dx(i,1) = dN_dxi(i,1)*JACOB_INV(1,1)+dN_dxi(i,2)*JACOB_INV(2,1)
              dN_dx(i,2) = dN_dxi(i,1)*JACOB_INV(1,2)+dN_dxi(i,2)*JACOB_INV(2,2)
          END DO

!        --- Assemble B-matrix ---
          BMAT_U(:,:) = 0.0D0
          DO i=1,NNODE
              BMAT_U(1,2*i-1) = dN_dx(i,1)
              BMAT_U(2,2*i  ) = dN_dx(i,2)
              BMAT_U(3,2*i-1) = dN_dx(i,2)
              BMAT_U(3,2*i  ) = dN_dx(i,1)
          END DO

!        --- Compute strain ---
          eps(:,1) = 0.0D0
          DO i=1,3
              DO j=1,2*NNODE
                  eps(i,1) = eps(i,1) + BMAT_U(i,j) * U(j)
              END DO
          END DO

!        --- Elastic stress ---
          CALL ELASTIC_STRESS(STRS=STRESS, DDSDDE=DDSDDE, EPS=eps(:,1), PROPS=PROPS)

!        --- Integration weight ---
          weight = gw(INTPT) * gw(INTPT) * JACOB_DET  ! simplified for 2D

!        --- Internal force vector ---
          DO i=1,2*NNODE
              DO j=1,3
                  Fint_u(i) = Fint_u(i) + BMAT_U(j,i)*STRESS(j)*weight
              END DO
          END DO

!        --- Element stiffness ---
          DO i=1,2*NNODE
              DO j=1,2*NNODE
                  DO k=1,3
                      DO l=1,3
                          AMATRX(i,j) = AMATRX(i,j) + BMAT_U(k,i)*DDSDDE(k,l)*BMAT_U(l,j)*weight
                      END DO
                  END DO
              END DO
          END DO

      END DO ! Gauss points

!     --- Assemble RHS ---
      DO i=1,2*NNODE
          RHS(i,1) = -Fint_u(i)
      END DO

      RETURN
END SUBROUTINE UEL

      !-------------------------------------------
! Helper: Plane stress / strain elasticity
SUBROUTINE ELASTIC_STRESS(STRS, DDSDDE, EPS, PROPS)
      DOUBLE PRECISION STRS(3), DDSDDE(3,3), EPS(3)
      DOUBLE PRECISION E, NU, D
      E  = PROPS(1)
      NU = PROPS(2)
      D  = 1.0D0
!     Plane strain
      DDSDDE(1,1) = E*(1.0D0-NU)/((1.0D0+NU)*(1.0D0-2.0D0*NU))
      DDSDDE(1,2) = E*NU/((1.0D0+NU)*(1.0D0-2.0D0*NU))
      DDSDDE(1,3) = 0.0D0
      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(2,2) = DDSDDE(1,1)
      DDSDDE(2,3) = 0.0D0
      DDSDDE(3,1) = 0.0D0
      DDSDDE(3,2) = 0.0D0
      DDSDDE(3,3) = E/(2.0D0*(1.0D0+NU))
      STRS(:) = MATMUL(DDSDDE, EPS)
      RETURN
END SUBROUTINE ELASTIC_STRESS


!-------------------------------------------
SUBROUTINE UEL_ERROR_MESSAGE(MSG)
      CHARACTER*(*) MSG
      WRITE(*,*) 'UEL ERROR:', MSG
      STOP
      RETURN
END SUBROUTINE

      SUBROUTINE SHAPE_FUNCTN(N, dN_dxi, H, dH_dxi, XI)

! Define input/output variables
! -----------------------------
    DIMENSION N(1,8), dN_dxi(8,2), H(1,4), dH_dxi(4,2)
    DIMENSION XI(2)

    DOUBLE PRECISION :: N, dN_dxi, H, dH_dxi
    DOUBLE PRECISION :: XI

! Shape function for 8 noded element (Q8)
! ---------------------------------------
    ! Serendipity element Q8 ---> N
    ! Corner nodes
    N = 0.0D0

    N(1,1) = -0.25D0 * (1.0D0 - XI(1)) * (1.0D0 - XI(2)) * (1.0D0 + XI(1) + XI(2))
    N(1,2) = -0.25D0 * (1.0D0 + XI(1)) * (1.0D0 - XI(2)) * (1.0D0 - XI(1) + XI(2))
    N(1,3) =  0.25D0 * (1.0D0 + XI(1)) * (1.0D0 + XI(2)) * (XI(1) + XI(2) - 1.0D0)
    N(1,4) = -0.25D0 * (1.0D0 - XI(1)) * (1.0D0 + XI(2)) * (1.0D0 + XI(1) - XI(2))
    ! Midside nodes
    N(1,5) = 0.5D0 * (1.0D0 - XI(1)**2) * (1.0D0 - XI(2))
    N(1,6) = 0.5D0 * (1.0D0 + XI(1)) * (1.0D0 - XI(2)**2)
    N(1,7) = 0.5D0 * (1.0D0 - XI(1)**2) * (1.0D0 + XI(2))
    N(1,8) = 0.5D0 * (1.0D0 - XI(1)) * (1.0D0 - XI(2)**2)

! Derivative of shape function, dN_dxi
! ------------------------------------
    dN_dxi = 0.0D0

    ! Derivative w.r.t XI(1) --> xi
    dN_dxi(1,1) = -0.25D0 * (-2.0D0 * XI(1) - XI(2)) * (1.0D0 - XI(2))
    dN_dxi(2,1) = -0.25D0 * (-2.0D0 * XI(1) + XI(2)) * (1.0D0 - XI(2))
    dN_dxi(3,1) =  0.25D0 *  (2.0D0 * XI(1) + XI(2)) * (1.0D0 + XI(2))
    dN_dxi(4,1) = -0.25D0 * (-2.0D0 * XI(1) + XI(2)) * (1.0D0 + XI(2))

    dN_dxi(5,1) = 0.5D0 * (-2.0D0 * XI(1)) * (1.0D0 - XI(2))
    dN_dxi(6,1) = 0.5D0 * (1.0D0 - XI(2)**2)
    dN_dxi(7,1) = 0.5D0 * (-2.0D0 * XI(1)) * (1.0D0 + XI(2))
    dN_dxi(8,1) = 0.5D0 * (-1.0D0 + XI(2)**2)

    ! Derivative w.r.t XI(2) --> eta
    dN_dxi(1,2) = -0.25D0 * (-XI(1) - 2.0D0 * XI(2)) * (1.0D0 - XI(1))
    dN_dxi(2,2) = -0.25D0 *  (XI(1) - 2.0D0 * XI(2)) * (1.0D0 + XI(1))
    dN_dxi(3,2) =  0.25D0 *  (XI(1) + 2.0D0 * XI(2)) * (1.0D0 + XI(1))
    dN_dxi(4,2) = -0.25D0 *  (XI(1) - 2.0D0 * XI(2)) * (1.0D0 - XI(1))

    dN_dxi(5,2) = 0.5D0 * (-1.0D0 + XI(1)**2)
    dN_dxi(6,2) = 0.5D0 * (-2.0D0 * XI(2)) * (1.0D0 + XI(1))
    dN_dxi(7,2) = 0.5D0 * (1.0D0 - XI(1)**2)
    dN_dxi(8,2) = 0.5D0 * (-2.0D0 * XI(2)) * (1.0D0 - XI(1))

! Shape function for 4 noded element (Q4)
! ---------------------------------------
    H = 0.0D0

    H(1,1) = 0.25D0 * (1.0D0 - XI(1)) * (1.0D0 - XI(2))
    H(1,2) = 0.25D0 * (1.0D0 + XI(1)) * (1.0D0 - XI(2))
    H(1,3) = 0.25D0 * (1.0D0 + XI(1)) * (1.0D0 + XI(2))
    H(1,4) = 0.25D0 * (1.0D0 - XI(1)) * (1.0D0 + XI(2))

! Derivative of shape function, dH_dxi
! ------------------------------------
    dH_dxi = 0.0D0

    ! Derivative w.r.t XI(1) --> xi
    dH_dxi(1,1) = -0.25D0 * (1.0D0 - XI(2))
    dH_dxi(2,1) =  0.25D0 * (1.0D0 - XI(2))
    dH_dxi(3,1) =  0.25D0 *: (1.0D0 + XI(2))
    dH_dxi(4,1) = -0.25D0 * (1.0D0 + XI(2))

    ! Derivative w.r.t XI(2) --> eta
    dH_dxi(1,2) = -0.25D0 * (1.0D0 - XI(1))
    dH_dxi(2,2) = -0.25D0 * (1.0D0 + XI(1))
    dH_dxi(3,2) =  0.25D0 * (1.0D0 + XI(1))
    dH_dxi(4,2) =  0.25D0 * (1.0D0 - XI(1))

    RETURN
END SUBROUTINE SHAPE_FUNCTN
