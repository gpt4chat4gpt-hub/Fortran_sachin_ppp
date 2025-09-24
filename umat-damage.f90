!DEC$ FREEFORM
! -------------------------------------------------------------------------------------------
! ABAQUS UMAT implementation for Transversely Isotropic Graphite flake with Gradient enhanced
! damage model
!
! UMAT file name  : umat-damage.f90
! Input file name :
! Date            : 20-07-2025
! -------------------------------------------------------------------------------------------
! Last edited on  : 18-09-2025
! Changes required: 
! Changes done    : Comments and description box of subroutine added.
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------  


SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
    RPL,DDSDDT,DRPLDE,DRPLDT, &
    STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
    NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
    CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
 
! ===========================================================================================
! This UMAT subroutine constructs the material constitutive law for transversely isotropic 
! graphite with the implementation of Gradient enhanced damage model.
! 
! Input to subroutine:-
! The material properties:
!   PROPS(1) ----> E11, Young's modulus of graphite in the out-of-plane direction
!   PROPS(2) ----> E22, Young's modulus of graphite in the in-plane direction
!   PROPS(3) ----> NU23, Poisson's ratio for in-plane direction
!   PROPS(4) ----> NU21, Poisson's ratio for out-of-plane direction
!   PROPS(5) ----> G21, Shear modulus
!   PROPS(6) ----> Damage threshold parameter (critical volumetric strain)
!   PROPS(7) ----> Minimum degradation function value
! 
! Output from subroutine (Return Variables):-
!   DDSDDE
!   STRESS
!   DDSDDT
!   RPL
!   DRPLDE
!   DRPLDT
! 
! State variables defined in the subroutine:-
!   STATEV(1) ----> Scalar damage variable
!   STATEV(2) ----> Non-local volumetric strain
! 
! ===========================================================================================

    INCLUDE 'ABA_PARAM.INC'      
    
    CHARACTER*80 CMNAME
    
! Define input/output variables
! -----------------------------
    INTEGER :: NDI, NSHR, NTENS, NSTATV, NPROPS
    INTEGER :: NOEL, NPT, LAYER, KSPT, JSTEP, KINC
    
    DOUBLE PRECISION :: STRESS, STATEV, DDSDDE, DDSDDT, DRPLDE
    DOUBLE PRECISION :: STRAN, DSTRAN, TIME, PREDEF, DPRED
    DOUBLE PRECISION :: PROPS, COORDS, DROT, DFGRD0, DFGRD1, PNEWDT
    DOUBLE PRECISION :: DTIME, TEMP, DTEMP, CELENT
    DOUBLE PRECISION :: SSE, SPD, SCD, RPL, DRPLDT
    
    DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
    DIMENSION DDSDDT(NTENS), DRPLDE(NTENS)
    DIMENSION STRAN(NTENS), DSTRAN(NTENS)
    DIMENSION TIME(2), PREDEF(1), DPRED(1)
    DIMENSION PROPS(NPROPS), COORDS(3), DROT(3,3)
    DIMENSION DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)

! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    DOUBLE PRECISION :: tolz = 1E-6
    
! Define internal variables
! -------------------------
    INTEGER :: i, j
    
    DOUBLE PRECISION :: E11, E22, E33, NU12, NU23, NU21, NU31, G21
    DOUBLE PRECISION :: gamma_
    DOUBLE PRECISION :: C(6,6), C11, C22, C33, C12, C13, C23, C44, C66
    
    DOUBLE PRECISION :: CMAT, I1_eps, dI1_eps, CE, kron_delta
    DOUBLE PRECISION :: dmin, eps_v_c, eps_v, epsbar_v_n1, epsbar_v_n
    DOUBLE PRECISION :: fd, dfd_dkappa, kappa, D, dD_dkappa, D_prev, D_new
    
    DIMENSION CMAT(NTENS,NTENS), CE(NTENS), kron_delta(NTENS)
    DIMENSION dI1_eps(NTENS), Sigma_n(NTENS)

! Assign Material properties (Direction: x (or 1) ---> out-of-plane; y & z (or 2  & 3) ---> in-plane)
! ---------------------------------------------------------------------------------------------------
    E11 =  PROPS(1)                   ! Young's modulus in the out-of-plane direction ('c' direction)
    E22 =  PROPS(2)                   ! Young's modulus in the in-plane direction ('a' direction)
    NU23 = PROPS(3)                   ! Poisson's ratio in in-plane direction 
    NU21 = PROPS(4)                   ! Poisson's ratio in out-of-plane direction
    G21 = PROPS(5)                    ! Shear modulus
    
    ! Degradation function properties
    eps_v_c = PROPS(6)                ! Damage threshold parameter (critical volumetric strain)
    dmin = PROPS(7)                   ! Minimum degradation function value
    
    E33 = E22                         ! Young's modulus in the in-plane direction ('a' direction)
    NU31 = NU21                       ! Poisson's ratio for out-of-plane direction
    NU12 = (NU21 * E11) / E33         ! Poisson's ratio for out-of-plane direction              
    
! Check for number of material parameters
! ---------------------------------------
    IF (NPROPS == 7) THEN
    ELSE
        PNEWDT 	= One_R
		IF (NOEL==1 .AND. NPT==1 .AND. KINC==1) THEN
			WRITE(6,*) "*******INPUTERROR !!! The number of material parameters are wrong!"
            		WRITE(6,*) "Check the material property input!"
			WRITE(6,*) "NUMBER ALLOWED: 7"
		END IF
    END IF

! Initialize variables
! --------------------
    CE(:) = Zero_R;         C(:,:) = Zero_R    
    DDSDDE(:,:) = Zero_R;   DDSDDT(:) = Zero_R
    DRPLDE(:) = Zero_R;     DRPLDT = Zero_R
    kron_delta(:) = Zero_R; dI1_eps(:) = Zero_R
    Sigma_n(:) = Zero_R
    
! Compute the stiffness matrix entries for transverly isotropic graphite flake 
! ----------------------------------------------------------------------------
    gamma_ = One_R / (One_R - NU23**2 - (Two_R * NU21 * NU12) - (Two_R * NU23 * NU21 * NU12))
    
    C11 = E11 * (One_R - NU23**2) * gamma_     
    C22 = E22 * (One_R - NU21 * NU12) * gamma_  
    C33 = C22                                   
    C23 = E22 * (NU23 + NU21 * NU12) * gamma_  
    C13 = E11 * (NU21 + NU23 * NU21) * gamma_   
    C12 = C13                                  
    C44 = G21                                  
    C66 = (C22 - C23) / Two_R               
    
! Formulate stiffness matrix
! --------------------------
    C(1,1) = C11
    C(1,2) = C12
    C(1,3) = C13

    C(2,1) = C12
    C(2,2) = C22
    C(2,3) = C23

    C(3,1) = C13
    C(3,2) = C23
    C(3,3) = C33

    C(4,4) = C44
    C(5,5) = C44
    C(6,6) = C66

! Declare state variables
! -----------------------
    ! Scalar damage variable
    D_prev = STATEV(1)
    
    ! Non-local volumetric strain
    epsbar_v_n = STATEV(2)
    
! Update total current strain
! ---------------------------
    DO i = 1, NTENS
    	STRAN(i) = STRAN(i) + DSTRAN(i)
    END DO
    
! Assign to tangent stiffness matrix
! ----------------------------------
    DO i = 1, NTENS                      
        DO j = 1, NTENS
            CMAT(i,j) = Zero_R
        END DO
    END DO
    
    ! Plane stress case (NDI==2 and NSHR==1)
    IF (NTENS == 3) THEN                 
    	DO i = 1, NDI
    		DO j = 1, NDI
    			CMAT(i,j) = C(i,j) - (C(i,3) * C(j,3)) / C(3,3)
    		END DO
    	END DO
    	
    	! Influence of shear modulus			
        i = NDI+1
        CMAT(i,i) = C(4,4)
        
    ! Plane strain case (NDI==3 and NSHR==1)
    ELSEIF (NTENS == 4) THEN
    	DO i = 1, NDI
    		DO j = 1, NDI
    			CMAT(i,j) = C(i,j)
    		END DO
    	END DO
    	
    	! Influence of shear modulus
    	i = NDI+1
    	CMAT(i,i) = C(4,4)
    	
    ELSE
    	! 3 dimensional case (NDI==3 and NSHR==3)
    	DO i = 1, NDI
        	DO j = 1, NDI
            		CMAT(i,j) = C(i,j)
        	END DO
    	END DO
	
	! Influence of shear modulus    
        DO i = NDI+1, NTENS                 
            CMAT(i,i) = C(i,i)
        END DO
      
    END IF

! Formulate kronecker delta
! -------------------------
    ! Plane stress case
    IF (NTENS <= 3) THEN
    	kron_delta(1) = One_R
        kron_delta(2) = One_R
        kron_delta(3) = Zero_R
        
    ! Plane strain case
    ELSEIF (NTENS == 4) THEN
    	kron_delta(1) = One_R
        kron_delta(2) = One_R
        kron_delta(3) = One_R
        kron_delta(4) = Zero_R
    
    ! 3 dimensional case
    ELSE 
    	kron_delta(1) = One_R
        kron_delta(2) = One_R
        kron_delta(3) = One_R
        kron_delta(4) = Zero_R
        kron_delta(5) = Zero_R
        kron_delta(6) = Zero_R
    END IF

! First invariant of strain tensor
! --------------------------------
    ! Plane stress case
    IF (NTENS <= 3) THEN
    	! First invariant of strain tensor
        I1_eps = (One_R - (CMAT(1,3) / CMAT(3,3))) * STRAN(1) + (One_R - (CMAT(2,3) / CMAT(3,3))) * STRAN(2)
        
        ! Derivative of first invariant of strain tensor w.r.t strain tensor
        DO i = 1, NTENS
        	dI1_eps(i) = (One_R - (CMAT(i,3) / CMAT(3,3))) * kron_delta(i)
        END DO
    
    ! Plane strain case
    ELSEIF (NTENS == 4) THEN
    	I1_eps = Zero_R
    	DO i = 1, NTENS
    		I1_eps = I1_eps + STRAN(i) * kron_delta(i)
    	END DO
   
    ! 3 dimensional case
    ELSE     
    	I1_eps = Zero_R
    	DO i = 1, NTENS
    		I1_eps = I1_eps + STRAN(i) * kron_delta(i)
    	END DO
    END IF

! Local damage variable / volumetric strain
! -----------------------------------------
    eps_v = I1_eps / Three_R
    
! Declare the non-local damage variable / non-local volumetric strain
! -------------------------------------------------------------------
    epsbar_v_n = TEMP
    epsbar_v_n1 = TEMP + DTEMP
    
    ! Define kappa
    kappa = epsbar_v_n1
  
! Call degradation function and compute the scalar damage variable
! ----------------------------------------------------------------
    CALL DEGRADATN_FUNCTN(dmin, eps_v_c, kappa, fd, dfd_dkappa)
    D_new = One_R - fd

! Ensure the irreversibility
! --------------------------          
    D = MAX(D_prev, D_new)
    
! Compute total current stress 
! ----------------------------
    ! Compute (CMAT : strain tensor)
    DO i = 1, NTENS
    	DO j = 1, NTENS
    		CE(i) = CE(i) + CMAT(i, j) * STRAN(j)
    	END DO
    END DO
    
    Sigma_n(:) = STRESS(:)
    ! Compute stress
    STRESS(:) = Zero_R 
    IF (ABS(D - D_prev) < 1.0e-5) THEN
        DO i = 1, NTENS
    	        STRESS(i) = Sigma_n(i)
        END DO
    ELSE
	    IF (NTENS <= 3) THEN
	    	DO i = 1, NTENS
	    		DO j = 1, NTENS
	    			STRESS(i) = STRESS(i) + (One_R - D) * CMAT(i, j) * STRAN(j) - (epsbar_v_n1 - eps_v) * (dI1_eps(i) / Three_R)
	    		END DO
	    	END DO
	    ELSE
	    	DO i = 1, NTENS
	    		DO j = 1, NTENS
	    			STRESS(i) = STRESS(i) + (One_R - D) * CMAT(i, j) * STRAN(j) - (epsbar_v_n1 - eps_v) * (kron_delta(i) / Three_R)
	    		END DO
	    	END DO
	    END IF
    END IF

! Assign tangent stiffness DDSDDE
! -------------------------------
    IF (NTENS <= 3) THEN
    	DO i = 1, NTENS
		DO j = 1, NTENS
	    		DDSDDE(i,j) = (One_R - D) * CMAT(i,j) + (dI1_eps(i) / Three_R) * (dI1_eps(j) / Three_R)
	    	END DO
	END DO
    ELSE
	DO i = 1, NTENS
		DO j = 1, NTENS
	    		DDSDDE(i,j) = (One_R - D) * CMAT(i,j) + (kron_delta(i) / Three_R) * (kron_delta(j) / Three_R)
	    	END DO
	END DO
    END IF	
    
! Assign source term RPL
! ----------------------
    RPL = eps_v - epsbar_v_n1
    
! Assign sensitivity DRPLDT
! -------------------------
    DRPLDT = - One_R
    
! Assign sensitivity DRPLDE
! -------------------------
    ! Plane stress case
    IF (NTENS <= 3) THEN
    	DO i = 1, NTENS
    		DRPLDE(i) = dI1_eps(i) / Three_R
    	END DO
    ELSE
    	DO i = 1, NTENS
    		DRPLDE(i) = kron_delta(i) / Three_R
    	END DO
    END IF
    
! Assign sensitivity DDSDDT
! -------------------------
    dD_dkappa = - dfd_dkappa
    
    IF (NTENS <= 3) THEN
    	DO i = 1, NTENS
    			DDSDDT(i) = -dD_dkappa * CE(i) - dI1_eps(i) / Three_R
	END DO
    ELSE
    	DO i = 1, NTENS
			DDSDDT(i) = -dD_dkappa * CE(i) - kron_delta(i) / Three_R
	END DO
    END IF

! Elastic strain energy
! ---------------------
    SSE = One_R/Two_R * DOT_PRODUCT(STRESS, DSTRAN)
    
! Store the STATEV 
! -----------------
    STATEV(1) = D
    STATEV(2) = kappa
    

! --------------------DEBUGGING PART------------------------
! Print results
! ----------------------------------------------------------
    IF(NOEL==1 .AND. NPT==1 .AND. KINC==1) THEN
    	!WRITE(6,*) '===== DDSDDE MATRIX ====='
        !DO i = 1, NTENS
        !    WRITE(6,'(A,I2,A,100(F12.2,1X))') 'DDSDDE(', i, ',:) = ', (DDSDDE(i,j), j=1,NTENS)
        !END DO
       ! 
     !	 WRITE(6,*) '===== STRESS VECTOR ====='
     !   DO i = 1, NTENS
     !       WRITE(6,'(A,I2,A,F12.2)') 'STRESS(', i, ') = ', STRESS(i)
     !   END DO
     !   
     !   WRITE(6,*) '===== CMAT MATRIX ====='
     !   DO i = 1, NTENS
     !       WRITE(6,'(A,I2,A,100(F12.2,1X))') 'CMAT(', i, ',:) = ', (CMAT(i,j), j=1,NTENS)
     !   END DO
    	
    END IF
    
    IF(NOEL==1 .AND. NPT==1) THEN
    	WRITE(6,*) ' '
    	WRITE(6,*) '===========================UMAT DEBUG OUTPUT=================================='
    	WRITE(6,*) 'KINC=', KINC, 'Kappa: ', kappa, 'D : ', D
    	!WRITE(6,*) 'STEP=', JSTEP
    	!WRITE(6,*) 'KINC=', KINC, ' ELEMENT=', NOEL, ' INTEGRATION PT=', NPT
    	!WRITE(6,*) ' ' 
    	!WRITE(6,*) 'TIME(1)=', TIME(1), ' TIME(2)=', TIME(2), ' DTIME=', DTIME
    	!WRITE(6,*) ' '
    	WRITE(6,*) '------------------- Local and non local Damage variables ---------------------'
    	WRITE(6,*) 'Threshold value eps_v_c=', eps_v_c
    	WRITE(6,*) 'Local damage variable or eps_v=', eps_v
    	WRITE(6,*) ' '
    	WRITE(6,*) 'Epsbar_v_n (prev)=', epsbar_v_n 
    	WRITE(6,*) 'TEMP or epsbar_v_n1=', epsbar_v_n1, ' DTEMP=', DTEMP
    	!WRITE(6,*) ' '
    	!write(6,*) 'kappa=', kappa
    	!WRITE(6,*) ' '
    	WRITE(6,*) '----------------------------- Damage functoin ---------------------------------'
    	WRITE(6,*) 'Degradation function fd = ', fd
    	!WRITE(6,*) 'Derivation of Degradation function dfd_dkappa = ', dfd_dkappa
    	!WRITE(6,*) 'dD_dkappa = ', dD_dkappa
    	!WRITE(6,*) ' '
    	! Strain and incremental strain
    	WRITE(6,*) '------------------- Strain ---------------------'
    	WRITE(6,*) 'Strain (STRAN):', (STRAN(i), i=1,NTENS)
    	WRITE(6,*) 'Incre Strain (DSTRAN):', (DSTRAN(i), i=1,NTENS)
    	WRITE(6,*) ' '
    	
        WRITE(6,*) '===== CMAT MATRIX ====='
        DO i = 1, NTENS
            WRITE(6,'(A,I2,A,100(F12.2,1X))') 'CMAT(', i, ',:) = ', (CMAT(i,j), j=1,NTENS)
        END DO
    	! Stress
    	WRITE(6,*) '------------------- Stress ---------------------'
    	DO i=1,NTENS
    		WRITE(6,'(A,I5,A,F15.6)') 'STRESS(', i, ') = ', STRESS(i)
    	END DO
    	WRITE(6,*) ' '
        ! Stress_previous
        WRITE(6,*) '------------------- Stress previous ---------------------'
        DO i=1,NTENS
                WRITE(6,'(A,I5,A,F15.6)') 'Sigma_n(', i, ') = ',  Sigma_n(i)
        END DO
        WRITE(6,*) ' '
    	! Tangent stiffness
    	WRITE(6,*) '------------------- DDSDDE matrix ---------------------'
    	DO i = 1, NTENS
            WRITE(6,'(A,I2,A,100(F12.2,1X))') 'DDSDDE(',i, ',:) = ', (DDSDDE(i,j), j=1,NTENS)
        END DO
        WRITE(6,*) ' '
    	! State variables (D and kappa)
    	WRITE(6,*) '------------------- STATEV ---------------------'
    	WRITE(6,*) 'D=', STATEV(1), ' kappa=', STATEV(2)
    	WRITE(6,*) ' '
    	WRITE(6,*) '===========================UMAT DEBUG OUTPUT=================================='
    END IF
     
    RETURN
END SUBROUTINE UMAT
  
    
! 
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
! 
! UMATHT subroutine
! 
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------


SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG, &
    STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED, &
    CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT, &
    NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

! ===========================================================================================
! This UMATHT subroutine constructs the balance equation of damage field & its sensitivities 
! for implementing the Gradient enhanced damage model.
! 
! Input to subroutine:-
! The material properties:
!   PROPS(1) ----> C_d, Damage regularization parameter ---> Conductivity
!   PROPS(2) ----> Cp, Small numerical value ---> Heat capacity
! 
! Output from subroutine (Return Variables):-
!   U
!   DUDT
!   DUDG
!   FLUX
!   DFDT
!   DFDG
! 
! ===========================================================================================
    INCLUDE 'ABA_PARAM.INC'

    CHARACTER*80 CMNAME
    
    DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD)
    DIMENSION DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD)
    DIMENSION TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
    
! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Define internal variables
! -------------------------
    INTEGER :: i
    
    DOUBLE PRECISION :: C_d, Cp, dU

! Define properties
! -----------------
    ! Damage regularization parameter ---> Conductivity
    C_d = PROPS(1)
    
    ! Heat capacity ---> small numerical value           
    Cp = PROPS(2)
    
! Assign sensitivity DUDT
! ----------------------- 
    DUDT = Cp
    
! Assign sensitivity DUDG
! -----------------------
    DUDG = Zero_R
    
! Assign sensitivity U, Internal energy
! -------------------------------------
    dU = DUDT * DTEMP
    U = dU
    
! Assign sensitivity FLUX, DFDG, DFDT
! -----------------------------------
    DO i = 1, NTGRD
    	FLUX(i) = - (C_d * C_d) * DTEMDX(i)
    	DFDG(i,i) = - (C_d * C_d) * One_R
    	DFDT(i) = Zero_R
    END DO

    RETURN
END SUBROUTINE UMATHT


! 
! -------------------------------------------------------------------------------------------
! Degradation function subroutine for graphite
! -------------------------------------------------------------------------------------------


SUBROUTINE DEGRADATN_FUNCTN(dmin, eps_v_c, kappa, fd, dfd_dkappa)

! ===========================================================================================
! This DEGRADATN_FUNCTN subroutine constructs the degradation function for the Graphite flake.
! 
! Output from subroutine (Return Variables):-
!   fd          ----> Degradation function value
!   dfd_dkappa  ----> First derivative of the degradation function
! 
! ===========================================================================================

! Define input variables
! ----------------------
    DOUBLE PRECISION :: dmin, eps_v_c, kappa
    
! Define output variables
! -----------------------
    DOUBLE PRECISION :: fd, dfd_dkappa
    
! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Degradation function for graphite
! ---------------------------------
    fd = ((One_R - dmin) / Two_R) * (One_R + TANH(-Four_R * kappa / eps_v_c)) + dmin
    
! Derviative of degradation function w.r.t history parameter (d_fd/d_kappa)
! -------------------------------------------------------------------------
    dfd_dkappa = -(Two_R * (One_R - dmin) / eps_v_c) * (One_R - TANH(-Four_R * kappa / eps_v_c)**2)
    
    RETURN
END SUBROUTINE DEGRADATN_FUNCTN

!
! ------------------------------------------------------------------------------------------- 
! ------------------------------------------------------------------------------------------- 
