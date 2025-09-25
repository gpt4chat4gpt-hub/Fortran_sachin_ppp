SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS, &
    PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME, &
    KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF, &
    LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

    INCLUDE 'ABA_PARAM.INC'

! Define input/output variables
! -----------------------------
    INTEGER :: MLVARX, NRHS, NDOFEL, NPROPS, NSVARS
    INTEGER :: MCRD, NNODE, JTYPE, KSTEP, KINC, JELEM, NDLOAD, MDLOAD, NPREDF, NJPROP
    
    DOUBLE PRECISION :: RHS, AMATRX, PROPS
    DOUBLE PRECISION :: SVARS, ENERGY, COORDS, U
    DOUBLE PRECISION :: DU, V, A, TIME, PARAMS
    DOUBLE PRECISION :: JDLTYP, ADLMAG, DDLMAG
    DOUBLE PRECISION :: PREDEF, JPROPS
    
    LOGICAL :: LFLAGS
    
    DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(NPROPS)
    DIMENSION SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL)
    DIMENSION DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(3)
    DIMENSION JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
    DIMENSION PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

! Float numbers / parameters
! --------------------------
    INTEGER, PARAMETER :: NGP = 9                ! Gauss quadrature 3x3
    INTEGER, PARAMETER :: NQ4 = 8                ! Number of nodes: 8 noded element
    INTEGER, PARAMETER :: NCOMP = 3              ! e11, e22, e12
    
    INTEGER, PARAMETER :: NSDV_PER_GP = 9        ! Number of SDV per node so NSVARS = NGP * NSDV_PER_GP
    INTEGER, PARAMETER :: PSTRAN_PSTRESS = 0     ! Switch for Plane Strain = 0 and Plane Stress = 1 
    
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
    LOGICAL, PARAMETER :: DEBUG_ON = .false.      ! Debugging
    INTEGER :: ierr                              ! Debugging
    CHARACTER(len=120) :: fname                  ! Debugging

! Define internal variables
! -------------------------
    INTEGER :: i, j, k, l, INTPT

    DOUBLE PRECISION :: XI, gp, gw
    DOUBLE PRECISION :: N, H, dN_dxi, dH_dxi
    DOUBLE PRECISION :: JACOB, JACOB_DET, JACOB_INV, CO_FACTOR, BMAT_U, BMAT_k
    DOUBLE PRECISION :: dN_dx, dH_dx
    
    DOUBLE PRECISION :: U_x, DU_x, epsbar_v, depsbar_v!, U_n1, epsbar_v_n1
    DOUBLE PRECISION :: eps_tr, deps, eps_n, epsbar_v_NL, depsbar_v_NL, epsbar_v_NL_n
    DOUBLE PRECISION :: SIGMA, STRESS_3, DDSDDE_3x3, DDSDDT_3, DRPLDE_3
    
    DOUBLE PRECISION :: weight, Fint_u, R_k
    DOUBLE PRECISION :: K_uu, K_uk, K_ku, K_kk
    DOUBLE PRECISION :: grad_k
    
    DOUBLE PRECISION :: X, Y 
    
    DOUBLE PRECISION :: E11, E22, NU23, NU21, G21, eps_v_c, dmin, C_d, cp, vec1, vec2
    
    DIMENSION XI(2), gp(2,NGP), gw(2,NGP)
    DIMENSION N(1,NNODE), H(1,NQ4), dN_dxi(NNODE,2), dH_dxi(NQ4,2)
    DIMENSION JACOB(2,2), JACOB_INV(2,2), CO_FACTOR(2,2), BMAT_U(NCOMP,2*NNODE), BMAT_k(2,NQ4)
    DIMENSION dN_dx(NNODE,2), dH_dx(NQ4,2)
    
    DIMENSION U_x(2*NNODE,1), DU_x(2*NNODE,1), epsbar_v(NQ4,1), depsbar_v(NQ4,1)!, U_n1(2*NNODE,1), epsbar_v_n1(NQ4,1)
    DIMENSION eps_tr(NCOMP,NGP), deps(NCOMP,NGP), eps_n(NCOMP,NGP), epsbar_v_NL(1,NGP), depsbar_v_NL(1,NGP), epsbar_v_NL_n(1,NGP)
    DIMENSION SIGMA(NCOMP), STRESS_3(NCOMP,NGP), DDSDDE_3x3(NCOMP,NCOMP), DDSDDT_3(NCOMP), DRPLDE_3(NCOMP)
    
    DIMENSION Fint_u(2*NNODE,1), R_k(NQ4,1)
    DIMENSION K_uu(2*NNODE,2*NNODE), K_uk(2*NNODE,NQ4), K_ku(NQ4,2*NNODE), K_kk(NQ4,NQ4)
    DIMENSION grad_k(2)
    
    DIMENSION X(NNODE), Y(NNODE)
    
    DIMENSION vec1(3), vec2(3)
     
! Variables defined for UMAT interface call in UEL
! ------------------------------------------------
    ! Switch definition for variables/parameters for Plane stress and plane strain case
    INTEGER, PARAMETER :: NDI = MERGE(3, 2, PSTRAN_PSTRESS == 0)
    INTEGER, PARAMETER :: NSHR = 1
    INTEGER, PARAMETER :: NTENS = MERGE(4, 3, PSTRAN_PSTRESS == 0)  
    INTEGER, PARAMETER :: NSTATV = 2
    
    INTEGER :: NOEL, NPT, LAYER, KSPT, JSTEP
    
    DOUBLE PRECISION :: STRESS, STATEV, STATEV_GP, DDSDDE, DDSDDT, DRPLDE
    DOUBLE PRECISION :: STRAN, DSTRAN
    DOUBLE PRECISION :: COORDS_GP, DROT, DFGRD0, DFGRD1, PNEWDT
    DOUBLE PRECISION :: TEMP, DTEMP, CELENT
    DOUBLE PRECISION :: SSE, SPD, SCD, RPL, DRPLDT
     
    CHARACTER*80 CMNAME
    
    DIMENSION STRESS(NTENS), STATEV(NSTATV, NGP), STATEV_GP(NSTATV), DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS)
    DIMENSION STRAN(NTENS), DSTRAN(NTENS)
    DIMENSION COORDS_GP(3), DROT(3,3)
    DIMENSION DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)

! Unpack material properties
! --------------------------
    E11 =  PROPS(1)             
    E22 =  PROPS(2)             
    NU23 = PROPS(3)             
    NU21 = PROPS(4)             
    G21 = PROPS(5)                
    eps_v_c = PROPS(6)
    dmin = PROPS(7)
    C_d = PROPS(8)
    cp = PROPS(9)
    vec1(1) = 1 !PROPS(10)         ! Define orientation of local CS Vector 1 
    vec1(2) = 0 !PROPS(11)
    vec1(3) = 0 !PROPS(12)
    
    vec2(1) = 0 !PROPS(13)         ! Define orientation of local CS Vector 2 
    vec2(2) = 1 !PROPS(14)
    vec2(3) = 0 !PROPS(15)
    
! Define nodal variables: Displacement, U_x and U_y; Non-local volumetric strain, epsbar_v
! ----------------------------------------------------------------------------------------
    U_x = Zero_R
    DU_x = Zero_R
    epsbar_v = Zero_R
    depsbar_v = Zero_R
    
    DO i = 1, NNODE
    	U_x(2*i-1,1) = U(3*i-2)
    	U_x(2*i,1) = U(3*i-1)
    	
    	DU_x(2*i-1,1) = DU(3*i-2,1)
    	DU_x(2*i,1) = DU(3*i-1,1)
    END DO
    
    DO i = 1, NQ4
    	epsbar_v(i,1) = U(3*i)
	depsbar_v(i,1) = DU(3*i,1)
    END DO
    
! Initialize variables (every element)
! ------------------------------------ 
    RHS(:,1) = Zero_R
    AMATRX(:,:) = Zero_R
    Fint_u(:,1) = Zero_R
    K_uu(:,:) = Zero_R
    K_uk(:,:) = Zero_R
    K_ku(:,:) = Zero_R
    K_kk(:,:) = Zero_R
    R_k(:,1) = Zero_R
    
    ! Initialize umat variables
    STRESS_3(:,:) = Zero_R; SIGMA = Zero_R
    DDSDDE_3x3(:,:) = Zero_R!; DDSDDE(:,:) = Zero_R
    DDSDDT_3 = Zero_R!; DDSDDT = Zero_R 
    DRPLDE_3 = Zero_R!; DRPLDE = Zero_R
    
! Call Gauss quadrature subroutine
! --------------------------------
    CALL GAUSS_PT_3x3(gp, gw)
        
! Assign coordinate values
! ------------------------
    DO i = 1, NNODE
    	X(i) = COORDS(1,i)
    	Y(i) = COORDS(2,i)
    END DO
    
! Begin loop over integration points
! ----------------------------------
    DO INTPT = 1, NGP
    	
    	! Assign local coordinates (XI(1) --> xi and XI(2) --> eta) at the integration point
    	XI(1) = gp(1, INTPT)
    	XI(2) = gp(2, INTPT)
    	
    	! Call shape functions and it's derivatives subroutine
    	CALL SHAPE_FUNCTN(XI, N, dN_dxi, H, dH_dxi)
    	
    	! Calculate the Jacobian
    	JACOB(:,:) = Zero_R; CO_FACTOR(:,:) = Zero_R
        JACOB_INV(:,:) = Zero_R
        JACOB_DET = Zero_R
        
    	DO i = 1, 2
    		DO j = 1, NNODE
    			! Derivative of coordinate X w.r.t xi and eta i.e. J(1,1) and J(1,2) 
    			JACOB(1,i) = JACOB(1,i) + dN_dxi(j,i) * X(j)
    		END DO
    	END DO
    	
    	DO i = 1, 2
    		DO j = 1, NNODE
    			! Derivative of coordinate Y w.r.t xi and eta i.e. J(2,1) and J(2,2) 
    			JACOB(2,i) = JACOB(2,i) + dN_dxi(j,i) * Y(j)
    		END DO
    	END DO
    	
    	! Calculate determinant of Jacobian
    	JACOB_DET = JACOB(1,1) * JACOB(2,2) - JACOB(1,2) * JACOB(2,1)
    	
    	! check if determinant is Zero
    	IF (JACOB_DET <= Zero_R) THEN 
            	WRITE(6,*) 'Determinant of Jacobian: ', JACOB_DET
    		CALL UEL_ERROR_MESSAGE('*** Determinant of Jacobian is 0 or negative ***') 
    	END IF
    	
    	! Calculate inverse of Jacobian
    	CO_FACTOR(1,1) = JACOB(2,2)
    	CO_FACTOR(2,2) = JACOB(1,1)
    	CO_FACTOR(1,2) = -JACOB(1,2)
    	CO_FACTOR(2,1) = -JACOB(2,1)
    	
        DO i = 1, 2
            DO j = 1, 2
    	        JACOB_INV(i,j) = CO_FACTOR(i,j) / JACOB_DET
            END DO
        END DO
    	
    	! Evaluate dN_dx for global coordinate --> Q8
        dN_dx(:,:) = Zero_R
    	DO i = 1, 2
    		DO j = 1, NNODE
    			! Derivative of shape function w.r.t x(i = 1) and w.r.t y(i = 2) 
    			dN_dx(j,i) = dN_dxi(j,1) * JACOB_INV(1,i) + dN_dxi(j,2) * JACOB_INV(2,i)
    		END DO
        END DO
    	
    	! Evaluate dH_dx for global coordinate --> Q4
        dH_dx(:,:) = Zero_R
    	DO i = 1, 2
    		DO j = 1, NQ4
    			! Derivative of shape function w.r.t x(i = 1) and w.r.t y(i = 2)
    			dH_dx(j,i) = dN_dxi(j,1) * JACOB_INV(1,i) + dN_dxi(j,2) * JACOB_INV(2,i)
    		END DO
    	END DO
    	
    	! Formulate B Matrix: BMAT_U(3 x 2*8) for displacement, U_x and U_y
    	DO i = 1, 3
    		DO j = 1, 2*NNODE
    			BMAT_U(i,j) = Zero_R
    		END DO
    	END DO
    	
    	DO i = 1, NNODE
    		BMAT_U(1,2*i-1) = dN_dx(i,1)
    		BMAT_U(2,2*i) = dN_dx(i,2)
    		BMAT_U(3,2*i) = dN_dx(i,1)
    		BMAT_U(3,2*i-1) = dN_dx(i,2)
    	END DO
    	
    	! Formulate B Matrix: BMAT_k(2 x 4) for Non-local volumetric strain, epsbar_v
    	DO i = 1, 2
    		DO j = 1, NQ4
    			BMAT_k(i,j) = Zero_R
    		END DO
    	END DO
    	
    	DO i = 1, NQ4
    		BMAT_k(1,i) = dH_dx(i,1)
    		BMAT_k(2,i) = dH_dx(i,2)
    	END DO
    	
    	! Kinematics: Compute strains at gauss points
    	DO i = 1, NCOMP
    		eps_tr(i,INTPT) = Zero_R; deps(i,INTPT) = Zero_R; eps_n(i,INTPT) = Zero_R
    		DO j = 1, 2*NNODE
    			eps_tr(i,INTPT) = eps_tr(i,INTPT) + BMAT_U(i,j) * U_x(j,1)
    			deps(i,INTPT) = deps(i,INTPT) + BMAT_U(i,j) * DU_x(j,1)
    		END DO
    		eps_n(i,INTPT) = eps_tr(i,INTPT) - deps(i,INTPT)
    	END DO
    	
    	! Compute the non-local volumetric strain at gauss points
    	epsbar_v_NL(1,INTPT) = Zero_R
        depsbar_v_NL(1,INTPT) = Zero_R
        epsbar_v_NL_n(1,INTPT) = Zero_R
        
    	DO i = 1, NQ4
    		epsbar_v_NL(1,INTPT) = epsbar_v_NL(1,INTPT) + N(1,i) * epsbar_v(i,1)
            	depsbar_v_NL(1,INTPT) = depsbar_v_NL(1,INTPT) + N(1,i) * depsbar_v(i,1)
        END DO
        epsbar_v_NL_n(1,INTPT) = epsbar_v_NL(1,INTPT) - depsbar_v_NL(1,INTPT)
    	
    	! Mapping variable to UMAT variables
    	! ----------------------------------
    	
    	! Strains for plane strain case
    	IF (NTENS == 4) THEN
    		STRAN(1) = eps_n(1,INTPT)
    		STRAN(2) = eps_n(2,INTPT) 
    		STRAN(3) = Zero_R
    		STRAN(4) = eps_n(3,INTPT) 
    		
    		DSTRAN(1) = deps(1,INTPT)
    		DSTRAN(2) = deps(2,INTPT)
    		DSTRAN(3) = Zero_R
    		DSTRAN(4) = deps(3,INTPT)
    	! Strains for plane stress case
    	ELSEIF (NTENS <= 3) THEN
    		STRAN(1) = eps_n(1,INTPT)
    		STRAN(2) = eps_n(2,INTPT)
    		STRAN(3) = eps_n(3,INTPT)

    		DSTRAN(1) = deps(1,INTPT)
    		DSTRAN(2) = deps(2,INTPT)
    		DSTRAN(3) = deps(3,INTPT)
    	! 3 dimensional case
    	ELSE
    		CALL UEL_ERROR_MESSAGE('*** UEL is not formulated to work for 3D case: Check the NTENS value in the parameters section and switch case (Refer line 49 and line 99 - 101).  ***')
    	END IF
    	
    	! Non-local volumetric strain
    	TEMP = epsbar_v_NL_n(1,INTPT)
    	DTEMP = depsbar_v_NL(1,INTPT)
    	
    	! Gauss point coordinates
    	COORDS_GP = Zero_R
    	DO i = 1, NNODE
    		COORDS_GP(1) = COORDS_GP(1) + N(1,i) * COORDS(1,i)
    		COORDS_GP(2) = COORDS_GP(2) + N(2,i) * COORDS(2,i)
    	END DO
    	COORDS_GP(3) = Zero_R	
    		
    	! Other mapping variables
    	CMNAME = 'Damage_UEL-UMAT'
    	CELENT = CELENT
    	NOEL = JELEM
    	NPT = INTPT
    	LAYER = 1
    	KSPT = 0
    	JSTEP = KSTEP
    	KINC = KINC
    	
    	! State variables
    	DO i = 1, NSTATV
    		STATEV(i,INTPT) = SVARS(NSDV_PER_GP*(INTPT-1)+i)
    		STATEV_GP(i) = STATEV(i,INTPT)
        END DO
    	
        DO i=1,3
          DO j=1,3
            DROT(i,j)=0.D0; DFGRD0(i,j)=0.D0; DFGRD1(i,j)=0.D0
          END DO
        END DO
        
    	! Call Damage_UMAT subroutine
    	! ---------------------------
	    CALL UMAT(STRESS,STATEV_GP,DDSDDE,SSE,SPD,SCD, &
    		RPL,DDSDDT,DRPLDE,DRPLDT, &
    		STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
    		NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS_GP,DROT,PNEWDT, &
    		CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

    	! Update the state variables
    	DO i = 1, NSTATV
            STATEV(i, INTPT) = STATEV_GP(i)
            SVARS(NSDV_PER_GP*(INTPT-1)+i) = STATEV(i, INTPT)
    	END DO
    	
    	! Update Stress and Material tangents
    	IF (NTENS == 4) THEN
            CALL CONVERT4TO3(STRESS, SIGMA)
            
            CALL CONVERT4x4TO3x3(DDSDDE, DDSDDE_3x3)
            
            CALL CONVERT4TO3(DDSDDT, DDSDDT_3)
            
            CALL CONVERT4TO3(DRPLDE, DRPLDE_3)
            
            DO i = 1, NCOMP
                STRESS_3(i,INTPT) = SIGMA(i)
            END DO
        
    	ELSEIF (NTENS <= 3) THEN
    		DO i = 1, NCOMP
    			STRESS_3(i,INTPT) = STRESS(i)
                	DDSDDT_3(i) = DDSDDT(i)
                	DRPLDE_3(i) = DRPLDE(i)
                
    			DO j = 1, NCOMP
    				DDSDDE_3x3(i,j) = DDSDDE(i,j)
    			END DO
    		END DO
    	ELSE
    		CALL UEL_ERROR_MESSAGE('*** UEL is not formulated to work for 3D case. ***')
    	END IF

    	! Formulate the element stiffness tangents and residuals
    	! ------------------------------------------------------
    	
    	! Integration weight
    	weight = gw(1,INTPT) * gw(2,INTPT) * JACOB_DET
    	
    	! Formulate ---> Fint_u
        ! Fint_u(:,1) = Fint_u(:,1) + MATMUL(TRANSPOSE(BMAT_U), STRESS_3(:,INTPT) ) * weight
    	DO i = 1, 2*NNODE
    		DO j = 1, NCOMP
    			Fint_u(i,1) = Fint_u(i,1) + BMAT_U(j,i) * STRESS_3(j,INTPT) * weight
    		END DO
        END DO

    	! Formulate ---> K_uu
        ! K_uu = K_uu + weight * MATMUL(TRANSPOSE(BMAT_U), MATMUL(DDSDDE_3x3,BMAT_U))
    	DO i = 1, 2*NNODE
    		DO j = 1, 2*NNODE
    			DO k = 1, NCOMP
    				DO l = 1, NCOMP
    					K_uu(i,j) = K_uu(i,j) + BMAT_U(k,i) * DDSDDE_3x3(k,l) * BMAT_U(l,j) * weight
    				END DO
    			END DO
    		END DO
        END DO

    	! Formulate ---> K_uk
    	DO i = 1, 2*NNODE
            DO j = 1, NQ4
                DO k = 1, NCOMP
                    K_uk(i,j) = K_uk(i,j) +  BMAT_U(k,i) * DDSDDT_3(k) * N(1,j) *  weight
                END DO
            END DO
        END DO
        
    	! Formulate ---> K_ku
    	DO i = 1, NQ4
            DO j = 1, 2*NNODE
                DO k = 1, NCOMP
    	            K_ku(i,j) = K_ku(i,j) - N(1,i) * DRPLDE_3(k) * BMAT_U(k,j) * weight
                END DO
            END DO
        END DO
    	
    	! Formulate ---> K_kk
    	DO i = 1, NQ4
            DO j = 1, NQ4
            	K_kk(i,j) = K_kk(i,j) - (N(1,i) * DRPLDT * N(1,j)) * weight
                DO k = 1, 2
    	            K_kk(i,j) = K_kk(i,j) + (BMAT_k(k,i) * C_d**2 * BMAT_k(k,j)) * weight
                END DO
            END DO
        END DO
    	
    	! Formulate the residual for Damage equation ---> R_k or F_epsbar
    	grad_k = Zero_R
        
        DO i = 1, NQ4
            DO j = 1, 2
                grad_k(j) = grad_k(j) + BMAT_k(j,i) * epsbar_v(i,1)
            END DO
        END DO
        
        DO i = 1, NQ4
            DO j = 1, 2
                R_k(i,1) = R_k(i,1) + ((BMAT_k(j,i) * C_d**2 * grad_k(j)) - N(1,i) * RPL) * weight
            END DO
        END DO
          
        ! Update the SVARS
        ! ----------------
        
        ! Strains
        SVARS(NSDV_PER_GP*(INTPT-1)+3) = eps_tr(1,INTPT)
        SVARS(NSDV_PER_GP*(INTPT-1)+4) = eps_tr(2,INTPT)
        SVARS(NSDV_PER_GP*(INTPT-1)+5) = eps_tr(3,INTPT)
        
        ! Non-local volumetric strain
        SVARS(NSDV_PER_GP*(INTPT-1)+6) = epsbar_v_NL(1,INTPT)
        
        ! Stress
        SVARS(NSDV_PER_GP*(INTPT-1)+7) = STRESS_3(1,INTPT)
        SVARS(NSDV_PER_GP*(INTPT-1)+8) = STRESS_3(2,INTPT)
        SVARS(NSDV_PER_GP*(INTPT-1)+9) = STRESS_3(3,INTPT)
        
        ! --------------------DEBUGGING PART------------------------
        ! Print results
        ! ----------------------------------------------------------
        
        IF (DEBUG_ON) THEN
		IF (KINC <= 3 .AND. JELEM == 1 .AND. INTPT ==1) THEN
			WRITE(6,*) '====================== UEL DEBUG OUTPUT : START =============================='
			write(6,'(/,79("-"))')
			write(6,'("UEL DEBUG  KSTEP=",I5,"  KINC=",I5,"   INTPT=",I3)') KSTEP, KINC, INTPT
			WRITE(6,*) ' '
		 	WRITE(6,*) 'TIME(1)=', TIME(1), ' TIME(2)=', TIME(2), ' DTIME=', DTIME
		 	WRITE(6,*) ' '
			! --- local coords at IP
			!write(6,'("XI = (",ES13.6,", ",ES13.6,")")') XI(1), XI(2)

			! --- nodal coordinates (geometry)
			!do j = 1, NNODE
			!  write(6,'("Node",I3,": X=",ES13.6,"  Y=",ES13.6)') j, X(j), Y(j)
			!end do

			! --- shape functions (if you want them)
		 	!write(6,'("N(j) (geom/disp shapes):")')
		    	!do j = 1, NNODE
		    	!  write(6,'("  N(",I2,") = ",ES13.6)') j, N(1,j)
		    	!end do

		    	! --- parent-space derivatives dN/dxi, dN/deta
		    	!write(6,'("dN_dxi (parent-space):")')
		    	!do j = 1, NNODE
		    	!  write(6,'("  dN_dxi(",I2,") = [",ES13.6,", ",ES13.6,"]")') j, dN_dxi(j,1), dN_dxi(j,2)
		    	!end do

		    	! --- Jacobian
		    	!write(6,'("JACOB = [",ES13.6,", ",ES13.6,"; ",ES13.6,", ",ES13.6,"]")') &
		    	!	JACOB(1,1), JACOB(1,2), JACOB(2,1), JACOB(2,2)
		    	!write(6,'("detJ  = ",ES13.6)') JACOB_DET
		    
		    	! --- inverse Jacobian
		    	!write(6,'("JACOB_INV = [",ES13.6,", ",ES13.6,"; ",ES13.6,", ",ES13.6,"]")') &
		    	!     JACOB_INV(1,1), JACOB_INV(1,2), JACOB_INV(2,1), JACOB_INV(2,2)

		    	! --- global derivatives dN/dx, dN/dy (AFTER you compute dN_dx)
		    	!write(6,'("dN_dx (global-space):")')
		    	!do j = 1, NNODE
		    	!  write(6,'("  dN_dx(",I2,") = [",ES13.6,", ",ES13.6,"]")') j, dN_dx(j,1), dN_dx(j,2)
		    	!end do
		    	write(6,'(79("-"),/)')

		END IF
    
		IF(JELEM==1 .AND. KINC == 1 .AND. INTPT ==1) THEN
	    	    	WRITE(6,*) ' '
		 	!WRITE(6,*) '--------------------------------- B matrix ----------------------------------'
	    	    	!write(6,*) 'BMAT_U:'
	    	    	!DO i = 1, 3
	    		!	write(6,'(10F10.4)') (BMAT_U(i,j), j=1,2*NNODE)
		 	!END DO
		 	!WRITE(6,*) ' '
	    	    	!write(6,*) 'BMAT_k:'
	    	    	!DO i = 1, 2
	    		!	write(6,'(10F10.4)') (BMAT_k(i,j), j=1,NQ4)
		 	!END DO
		 	!WRITE(6,*) ' '
		    	WRITE(6,*) '---------------------------------- Strains ---------------------------------'
	    	    	write(6,*) 'Epsilon at GP ', INTPT, ':', (eps(i,INTPT), i=1,NCOMP)
	    	    	write(6,*) 'Depsilon at GP ', INTPT, ':', (deps(i,INTPT), i=1,NCOMP)
		    	WRITE(6,*) ' '
	    	    	write(6,*) 'epsbar_v_NL at GP ', INTPT, ':', epsbar_v_NL(1,INTPT)
		    	WRITE(6,*) ' '
		    	WRITE(6,*) '------------------------- Strains passed to umat ---------------------------'
	    	    	write(6,*) 'STRAN:', (STRAN(i), i=1,NTENS)
	    	    	write(6,*) 'DSTRAN:', (DSTRAN(i), i=1,NTENS)
		    	WRITE(6,*) ' '
		    	write(6,*) 'Stress computed by umat'
	    	    	write(6,*) 'STRESS:', (STRESS(i),i=1,NTENS)
	    	    	write(6,*) 'STATEV:', (STATEV_GP(i),i=1,NSTATV)
		    	WRITE(6,*) ' '
		    	write(6,*) 'Stress after umat call - resized'
	    	    	write(6,*) 'STRESS_3', INTPT, ':', (STRESS_3(i,INTPT),i=1,NCOMP)
		    	WRITE(6,*) ' '
	    	    	WRITE(6,*) '--------------------------- DDSDDE matrix -----------------------------------'
		    	write(6,*) 'DDSDDE computed by umat'
	    	    	DO i = 1, NTENS
		        	WRITE(6,'(A,I2,A,100(F12.2,1X))') 'DDSDDE(',i, ',:) = ', (DDSDDE(i,j), j=1,NTENS)
		    	END DO
		    	WRITE(6,*) ' '
		    	write(6,*) 'DDSDDE after umat call - resized'
		    	DO i = 1, NCOMP
		        	WRITE(6,'(A,I2,A,100(F12.2,1X))') 'DDSDDE_3x3(',i, ',:) = ', (DDSDDE_3x3(i,j), j=1,NCOMP)
		    	END DO
		    	WRITE(6,*) ' '
		    
		END IF
        END IF
    	
    END DO   ! End loop over integration points
    	
! Assembly into AMATRX and RHS
! ----------------------------
    
    ! Assemble K_uu
    DO i = 1, NNODE
        DO j = 1, NNODE
            AMATRX(3*i-2,3*j-2) = K_uu(2*i-1,2*j-1)
            AMATRX(3*i-1,3*j-2) = K_uu(2*i,2*j-1)
            AMATRX(3*i-2,3*j-1) = K_uu(2*i-1,2*j)
            AMATRX(3*i-1,3*j-1) = K_uu(2*i,2*j)
        END DO
    END DO

    ! Assemble K_uk
    DO i = 1, NNODE
        DO j = 1, NQ4
            AMATRX(3*i-2,3*j) = K_uk(2*i-1,j)
            AMATRX(3*i-1,3*j) = K_uk(2*i,j)
        END DO
    END DO

    ! Assemble K_ku
    DO i = 1, NQ4
        DO j = 1, NNODE
            AMATRX(3*i,3*j-2) = K_ku(i,2*j-1)
            AMATRX(3*i,3*j-1) = K_ku(i,2*j)
        END DO
    END DO

    ! Assemble K_kk
    DO i = 1, NQ4
        DO j = 1, NQ4
            AMATRX(3*i,3*j) = K_kk(i,j)
        END DO
    END DO
    
    ! Assemble RHS
    DO i = 1, NNODE
        !DO j = 1, NNODE
            RHS(3*i-2,1) = RHS(3*i-2,1) - Fint_u(2*i-1,1) !K_uu(2*i-1,j)*U_x(j,1)
            RHS(3*i-1,1) = RHS(3*i-1,1) - Fint_u(2*i,1) !K_uu(2*i,j)*U_x(j,1)
        !END DO
    END DO
    
    DO i = 1, NQ4
        RHS(3*i,1) = RHS(3*i,1) - R_k(i,1)
    END DO
     
       
! Debugging block for UEL: AMATRX, R, DOF mapping
! -----------------------------------------------
   IF (DEBUG_ON) THEN
	fname = './AMATRX_debug.dat'

	! Open file for writing
	OPEN(UNIT=20, FILE=fname, STATUS='REPLACE', IOSTAT=ierr)
	IF (ierr /= 0) THEN
     		WRITE(*,*) 'ERROR: Could not open debug file ', fname
	ELSE
	
     	! Header
		WRITE(20,*) '==============================='
		WRITE(20,*) 'UEL Debug - Increment=', KINC
		WRITE(20,*) '==============================='

     	! Print DOF mapping
     		WRITE(20,*) 'Node DOF mapping:'
     		DO i = 1, NNODE
         		WRITE(20,'(A,I3,A,I3,A,I3)') 'Node ', i, ': Ux->', 3*i-2, ', Uy->', 3*i-1, ', epsbar_v->', 3*i
     		END DO

     	! Print AMATRX (stiffness) matrix
     		WRITE(20,*) ''
     		WRITE(20,*) 'AMATRX matrix:'
     		DO i = 1, NDOFEL
         		WRITE(20,'(100(E15.6,1X))') (AMATRX(i,j), j=1,NDOFEL)
     		END DO

     	! Print Residual vector
     		WRITE(20,*) ''
     		WRITE(20,*) 'Residual vector R:'
     		DO i = 1, NDOFEL
         		WRITE(20,'(E15.6)') RHS(i,1)
     		END DO

     	! Flush buffer
     	FLUSH(20)
     	CLOSE(20)
	END IF
	WRITE(6,*) '======================== UEL DEBUG OUTPUT : END =============================='
    END IF
	
    RETURN
END SUBROUTINE UEL

      !-------------------------------------------
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
    RPL,DDSDDT,DRPLDE,DRPLDT, &
    STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
    NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
    CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
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
    DOUBLE PRECISION :: vec1, vec2, ROTMTX_Q, C_global, QTQ
    
    DOUBLE PRECISION :: CMAT, I1_eps, dI1_eps, CE, kron_delta
    DOUBLE PRECISION :: dmin, eps_v_c, eps_v, epsbar_v_n1, epsbar_v_n
    DOUBLE PRECISION :: fd, dfd_dkappa, kappa, D, dD_dkappa, D_prev, D_new
    DOUBLE PRECISION :: Sigma_n, eps_tr
    
    DIMENSION vec1(3), vec2(3), ROTMTX_Q(3,3), C_global(6,6), QTQ(3,3)
    DIMENSION CMAT(NTENS,NTENS), CE(NTENS), kron_delta(NTENS)
    DIMENSION dI1_eps(NTENS), Sigma_n(NTENS), eps_tr(NTENS)
   
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

    ! Assign rotation vector
    vec1(1) = 1 !PROPS(10)            ! Define orientation of local CS Vector 1 
    vec1(2) = 0 !PROPS(11)
    vec1(3) = 0 !PROPS(12)
    
    vec2(1) = 0 !PROPS(13)            ! Define orientation of local CS Vector 2 
    vec2(2) = 1 !PROPS(14)
    vec2(3) = 0 !PROPS(15)
    
    E33 = E22                         ! Young's modulus in the in-plane direction ('a' direction)
    NU31 = NU21                       ! Poisson's ratio for out-of-plane direction
    NU12 = (NU21 * E11) / E33         ! Poisson's ratio for out-of-plane direction             

! Check for number of material parameters
! ---------------------------------------
    IF (NPROPS == 9) THEN
    ELSE
        PNEWDT 	= One_R
		IF (NOEL==1 .AND. NPT==1 .AND. KINC==1) THEN
			WRITE(6,*) "*******INPUTERROR !!! The number of material parameters are wrong!"
            		WRITE(6,*) "Check the material property input!"
			WRITE(6,*) "NUMBER ALLOWED: 9"
		END IF
    END IF
    
! Initialize variables
! --------------------
    CE(:) = Zero_R;         C(:,:) = Zero_R    
    DDSDDE(:,:) = Zero_R;   DDSDDT(:) = Zero_R
    DRPLDE(:) = Zero_R;     DRPLDT = Zero_R
    kron_delta(:) = Zero_R; dI1_eps(:) = Zero_R
    Sigma_n(:) = Zero_R
    C_global = Zero_R;      ROTMTX_Q = Zero_R
    
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
    
! Update total current strain ---> Global strains passed from UEL
! ---------------------------------------------------------------
    DO i = 1, NTENS
    	eps_tr(i) = STRAN(i) !+ DSTRAN(i)
    END DO

! Convert the Orientation from Local Coordinates to Global Coordinates
! --------------------------------------------------------------------
    ! Formulate the Rotation matrix
    CALL ROTATN_MATRIX(Vec1, Vec2, ROTMTX_Q)
    
    ! Perform the rotation of 4th order Tensor from local to global CS
    CALL TENSOR6x6_LOC_2_GLOBL(C, ROTMTX_Q, C_global)
    
    ! Compute Q^T * Q
    !do i=1,3
    !   do j=1,3
    !      QTQ(i,j) = ROTMTX_Q(1,i)*ROTMTX_Q(1,j) + ROTMTX_Q(2,i)*ROTMTX_Q(2,j) + ROTMTX_Q(3,i)*ROTMTX_Q(3,j)
    !   end do
    !end do
    !
    !! Print QTQ to check
    !write(*,*) 'Q^T Q = '
    !do i=1,3
    !   write(*,'(3F12.6)') QTQ(i,1), QTQ(i,2), QTQ(i,3)
    !end do
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
    			CMAT(i,j) = C_global(i,j) - (C_global(i,3) * C_global(j,3)) / C_global(3,3)
    		END DO
    	END DO
    	
    	! Influence of shear modulus			
        i = NDI+1
        CMAT(i,i) = C_global(4,4)
        
    ! Plane strain case (NDI==3 and NSHR==1)
    ELSEIF (NTENS == 4) THEN
    	DO i = 1, NDI
    		DO j = 1, NDI
    			CMAT(i,j) = C_global(i,j)
    		END DO
    	END DO
    	
    	! Influence of shear modulus
    	i = NDI+1
    	CMAT(i,i) = C_global(4,4)
    	
    ELSE
    	! 3 dimensional case (NDI==3 and NSHR==3)
    	DO i = 1, NDI
        	DO j = 1, NDI
            		CMAT(i,j) = C_global(i,j)
        	END DO
    	END DO
	
	! Influence of shear modulus    
        DO i = NDI+1, NTENS                 
            CMAT(i,i) = C_global(i,i)
        END DO
      
    END IF

! Formulate kronecker delta
! -------------------------
    ! Plane stress case
    IF (NTENS == 3) THEN
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
    IF (NTENS == 3) THEN
    	! First invariant of strain tensor
        I1_eps = (One_R - (CMAT(1,3) / CMAT(3,3))) * eps_tr(1) + (One_R - (CMAT(2,3) / CMAT(3,3))) * eps_tr(2)
        
        ! Derivative of first invariant of strain tensor w.r.t strain tensor
        DO i = 1, NTENS
        	dI1_eps(i) = (One_R - (CMAT(i,3) / CMAT(3,3))) * kron_delta(i)
        END DO
    
    ! Plane strain case
    ELSEIF (NTENS == 4) THEN
    	I1_eps = Zero_R
    	DO i = 1, NTENS
    		I1_eps = I1_eps + eps_tr(i) * kron_delta(i)
    	END DO
   
    ! 3 dimensional case
    ELSE     
    	I1_eps = Zero_R
    	DO i = 1, NTENS
    		I1_eps = I1_eps + eps_tr(i) * kron_delta(i)
    	END DO
    END IF

! Local damage variable
! ---------------------
    eps_v = I1_eps / Three_R
    
! Declare the non-local damage variable
! -------------------------------------
    epsbar_v_n = TEMP
    epsbar_v_n1 = TEMP !+ DTEMP
    
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
    		CE(i) = CE(i) + CMAT(i, j) * eps_tr(j)
    	END DO
    END DO
    
    Sigma_n(:) = STRESS(:)
    ! Compute stress
    STRESS(:) = Zero_R 
    IF (ABS(D - D_prev) < 1.0e-7) THEN
        DO i = 1, NTENS
    	        STRESS(i) = Sigma_n(i)
        END DO
    ELSE
	    IF (NTENS <= 3) THEN
	    	DO i = 1, NTENS
	    		DO j = 1, NTENS
	    			STRESS(i) = STRESS(i) + (One_R - D) * CMAT(i, j) * eps_tr(j) - (epsbar_v_n1 - eps_v) * (dI1_eps(i) / Three_R)
	    		END DO
	    	END DO
	    ELSE
	    	DO i = 1, NTENS
	    		DO j = 1, NTENS
	    			STRESS(i) = STRESS(i) + (One_R - D) * CMAT(i, j) * eps_tr(j) - (epsbar_v_n1 - eps_v) * (kron_delta(i) / Three_R)
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
	    		!IF (kappa >= eps_v_c .AND. kappa == epsbar_v_n1) THEN
	    		!	DDSDDE(i,j) = DDSDDE(i,j) + dfd_dkappa * CE(i) * (dI1_eps(j) / Three_R)
	    		!END IF
	    	END DO
	END DO
    ELSE
	DO i = 1, NTENS
		DO j = 1, NTENS
	    		DDSDDE(i,j) = (One_R - D) * CMAT(i,j) + (kron_delta(i) / Three_R) * (kron_delta(j) / Three_R)
	    		!IF (kappa >= eps_v_c .AND. kappa == epsbar_v_n1) THEN
	    		!	DDSDDE(i,j) = DDSDDE(i,j) + dfd_dkappa * CE(i) * (kron_delta(j) / Three_R)
	    		!END IF
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
    IF (NTENS == 3) THEN
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
    
    IF (NTENS == 3) THEN
    	DO i = 1, NTENS
    			DDSDDT(i) = -dD_dkappa * CE(i) - dI1_eps(i) / Three_R
	    END DO
    ELSE
    	DO i = 1, NTENS
			DDSDDT(i) = -dD_dkappa * CE(i) - kron_delta(i) / Three_R
	    END DO
    END IF
    
! Store the STATEV 
! -----------------
    STATEV(1) = D
    STATEV(2) = kappa
    
! --------------------DEBUGGING PART------------------------
! Print results
! ----------------------------------------------------------
    
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
    	WRITE(6,*) 'Strain (STRAN):', (eps_tr(i), i=1,NTENS)
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
       ! DO i=1,NTENS
    !		WRITE(6,'(A,I5,A,F15.6)') 'Sigma_n(', i, ') = ', Sigma_n(i)
    !	END DO
     !   WRITE(6,*) ' '
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

SUBROUTINE DEGRADATN_FUNCTN(dmin, eps_v_c, kappa, fd, dfd_dkappa)
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

SUBROUTINE GAUSS_PT_3x3(gp, gw)

! Define input/output variables
! -----------------------------
    DOUBLE PRECISION :: gp, gw, alpha
    
    DIMENSION gp(2,9), gw(2,9)

! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Define gauss quadrature points
! ------------------------------
    alpha = SQRT(Three_R/Five_R)
    
    gp(1,1) = - alpha       
    gp(1,2) = Zero_R
    gp(1,3) = alpha
    gp(1,4) = - alpha
    gp(1,5) = Zero_R
    gp(1,6) = alpha
    gp(1,7) = - alpha
    gp(1,8) = Zero_R
    gp(1,9) = alpha
    
    gp(2,1) = - alpha
    gp(2,2) = - alpha
    gp(2,3) = - alpha
    gp(2,4) = Zero_R
    gp(2,5) = Zero_R
    gp(2,6) = Zero_R
    gp(2,7) = alpha
    gp(2,8) = alpha
    gp(2,9) = alpha

! Define gauss quadrature weights
! -------------------------------
    gw(1,1) = Five_R/Nine_R
    gw(1,2) = Eight_R/Nine_R
    gw(1,3) = Five_R/Nine_R
    gw(1,4) = Five_R/Nine_R
    gw(1,5) = Eight_R/Nine_R
    gw(1,6) = Five_R/Nine_R
    gw(1,7) = Five_R/Nine_R
    gw(1,8) = Eight_R/Nine_R
    gw(1,9) = Five_R/Nine_R
    
    gw(2,1) = Five_R/Nine_R
    gw(2,2) = Five_R/Nine_R
    gw(2,3) = Five_R/Nine_R
    gw(2,4) = Eight_R/Nine_R
    gw(2,5) = Eight_R/Nine_R
    gw(2,6) = Eight_R/Nine_R
    gw(2,7) = Five_R/Nine_R
    gw(2,8) = Five_R/Nine_R
    gw(2,9) = Five_R/Nine_R
    
    RETURN
END SUBROUTINE GAUSS_PT_3x3

    
!
! -------------------------------------------------------------------------------------------
! Convert Matrix 4x4 to 3x3 matrix subroutine
! -------------------------------------------------------------------------------------------


SUBROUTINE CONVERT4x4TO3x3(MATRIX, MATRIX_3x3)

! Define input/output variables
! -----------------------------
    DOUBLE PRECISION :: MATRIX, MATRIX_3x3
    
    DIMENSION MATRIX(4,4), MATRIX_3x3(3,3)
    
! Assign values
! -------------
    MATRIX_3x3 = 0.0D0
    
    MATRIX_3x3(1,1) = MATRIX(1,1)
    MATRIX_3x3(1,2) = MATRIX(1,2)
    MATRIX_3x3(1,3) = MATRIX(1,4)
    
    MATRIX_3x3(2,1) = MATRIX(2,1)
    MATRIX_3x3(2,2) = MATRIX(2,2)
    MATRIX_3x3(2,3) = MATRIX(2,4)
    
    MATRIX_3x3(3,1) = MATRIX(4,1)
    MATRIX_3x3(3,2) = MATRIX(4,2)
    MATRIX_3x3(3,3) = MATRIX(4,4)

    RETURN
END SUBROUTINE CONVERT4x4TO3x3
  
    
!    
! -------------------------------------------------------------------------------------------
! Convert Vector(4) to Vector(3) subroutine
! -------------------------------------------------------------------------------------------

    
SUBROUTINE CONVERT4TO3(VECTOR, VECTOR_3)

! Define input/output variables
! -----------------------------
    DOUBLE PRECISION :: VECTOR, VECTOR_3
    
    DIMENSION VECTOR(4), VECTOR_3(3)
    
! Assign values
! -------------
    VECTOR_3 = 0.0D0
    
    VECTOR_3(1) = VECTOR(1)                
    VECTOR_3(2) = VECTOR(2)                
    VECTOR_3(3) = VECTOR(4)                

    RETURN
END SUBROUTINE CONVERT4TO3


! ---------------------------------------------------------------
! Rotation of 4th order Tensor from local to global CS subroutine
! ---------------------------------------------------------------


SUBROUTINE TENSOR6x6_LOC_2_GLOBL(C_local_VGT, ROTMTX_Q, C_global_VGT)

! Define input/output variables
! -----------------------------
    INTEGER :: i, j, k, l, p, q, r, s
    
    DOUBLE PRECISION :: C_local_VGT, ROTMTX_Q, C_global_VGT
    DOUBLE PRECISION :: C_local_TNSR, C_global_TNSR
    
    DIMENSION C_local_VGT(6,6), ROTMTX_Q(3,3), C_global_VGT(6,6)
    DIMENSION C_local_TNSR(3,3,3,3), C_global_TNSR(3,3,3,3)
    
! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Initialize Variables
! --------------------
    DO i = 1, 3
    	DO j = 1, 3
    		DO k = 1, 3
    			DO l = 1, 3
    				C_local_TNSR(i,j,k,l) = Zero_R
    				C_global_TNSR(i,j,k,l) = Zero_R
    			END DO
    		END DO
    	END DO
    END DO
    
! Transform local stiffness tensor from voigt to tensor
! -----------------------------------------------------
    CALL VOIGT6x6_2_TENSOR3x3x3x3(C_local_VGT, C_local_TNSR)
    
! Perform rotation of 4th order Tensor
! ------------------------------------
    DO i = 1, 3
    	DO j = 1, 3
    		DO k = 1, 3
    			DO l = 1, 3
    				DO p = 1, 3
    					DO q = 1, 3
    						DO r = 1, 3
    							DO s = 1, 3
    								C_global_TNSR(i,j,k,l) = C_global_TNSR(i,j,k,l) + ROTMTX_Q(i,p) * ROTMTX_Q(j,q) * ROTMTX_Q(k,r) * ROTMTX_Q(l,s) * C_local_TNSR(p,q,r,s)
    							END DO
    						END DO
    					END DO
    				END DO
    			END DO
    		END DO
    	END DO
    END DO
    
! Transform the rotated stiffness tensor (local --> global) from tensor to voigt
! ------------------------------------------------------------------------------
    CALL TENSOR3x3x3x3_2_VOIGT6x6(C_global_TNSR, C_global_VGT)
    
    RETURN
END SUBROUTINE TENSOR6x6_LOC_2_GLOBL

    
!
! -------------------------------------------------------------------------------------------
! Formulation of rotation matrix subroutine
! -------------------------------------------------------------------------------------------


SUBROUTINE ROTATN_MATRIX(Vector1, Vector2, ROTMTX_Q)
    
! Define input/output variables
! -----------------------------
    INTEGER :: i
        
    DOUBLE PRECISION :: Vector1, Vector2, ROTMTX_Q
    DOUBLE PRECISION :: Vec1, Vec2
    DOUBLE PRECISION :: Ut_vec1, Ut_vec2, Ut_vec3
    DOUBLE PRECISION :: Norm_1, Inner_prodt, Norm_2
    
    DIMENSION Vector1(3), Vector2(3), ROTMTX_Q(3,3)
    DIMENSION Vec1(3), Vec2(3)
    DIMENSION Ut_vec1(3), Ut_vec2(3), Ut_vec3(3)
    
! Normalize the Vector 1
! ----------------------
    Vec1 = Vector1
    
    ! Compute norm of vector 1
    CALL NORM_OF_VECT(Vec1, Norm_1)
    
    ! Normalize vector 1 to unit vector
    DO i = 1, 3
    	Ut_vec1(i) = Vec1(i) / Norm_1
    END DO
    
! Use provided vector 2, project it perpendicular to vec1 and normalize vector 2
! ------------------------------------------------------------------------------
    ! Project Vector 2 onto vector 1
    Inner_prodt = Vec1(1)*Vector2(1) + Vec1(2)*Vector2(2) + Vec1(3)*Vector2(3)
    
    ! Project it perpendicular to vector 1
    DO i = 1, 3
    	Vec2(i) = Vector2(i) - Inner_prodt * Ut_vec1(i)
    END DO
    
    ! Compute norm of vector 2
    CALL NORM_OF_VECT(Vec2, Norm_2)
    
    ! Normalize vector 2 to unit vector
    DO i = 1, 3
    	Ut_vec2(i) = Vec2(i) / Norm_2
    END DO
    
! Compute the vector 3 orthonormal to vector 1 and vector 3
! ---------------------------------------------------------
    ! Cross product of unit vectors 1 and 2
    Ut_vec3(1) = Ut_vec1(2)*Ut_vec2(3) - Ut_vec1(3)*Ut_vec2(2)
    Ut_vec3(2) = Ut_vec1(3)*Ut_vec2(1) - Ut_vec1(1)*Ut_vec2(3)
    Ut_vec3(3) = Ut_vec1(1)*Ut_vec2(2) - Ut_vec1(2)*Ut_vec2(1)
    
! Formulate the Rotation matrix
! -----------------------------
    ROTMTX_Q(:,1) = Ut_vec1
    ROTMTX_Q(:,2) = Ut_vec2
    ROTMTX_Q(:,3) = Ut_vec3
    
    RETURN
END SUBROUTINE ROTATN_MATRIX

    
!
! -------------------------------------------------------------------------------------------
! Norm of vector subroutine
! -------------------------------------------------------------------------------------------


SUBROUTINE NORM_OF_VECT(Vector, Norm_)

! Define input/output variables
! -----------------------------
    DOUBLE PRECISION :: Vector, Norm_
    
    DIMENSION Vector(3)
    
! Compute norm of a vector
! ------------------------
    Norm_ = SQRT(Vector(1)*Vector(1) + Vector(2)*Vector(2) + Vector(3)*Vector(3))
    
    RETURN
END SUBROUTINE NORM_OF_VECT

    
!
! -------------------------------------------------------------------------------------------
! Voigt to Tensor transformation subroutine
! -------------------------------------------------------------------------------------------


SUBROUTINE VOIGT6x6_2_TENSOR3x3x3x3(C_voigt, C_tensor)

! Define input/output variables
! -----------------------------
    INTEGER :: i, j, k, l, m, n
    INTEGER :: V2T
    
    DOUBLE PRECISION :: C_voigt, C_tensor
    DOUBLE PRECISION :: factor
    
    DIMENSION C_voigt(6,6), C_tensor(3,3,3,3)
    DIMENSION factor(6), V2T(6,2)

! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Assign values
! -------------
    ! For converting Engineering shear strain to true strain tensor
    factor(1) = One_R; factor(4) = Two_R
    factor(2) = One_R; factor(5) = Two_R
    factor(3) = One_R; factor(6) = Two_R
    
    ! Indices in Abaqus voigt notation
    V2T(1,1) = 1; V2T(1,2) = 1
    V2T(2,1) = 2; V2T(2,2) = 2
    V2T(3,1) = 3; V2T(3,2) = 3
    V2T(4,1) = 1; V2T(4,2) = 2
    V2T(5,1) = 1; V2T(5,2) = 3
    V2T(6,1) = 2; V2T(6,2) = 3
    
    ! Formulate the tensor C_tensor
    DO m = 1, 6
    	i = V2T(m,1)
    	j = V2T(m,2)
    	DO n = 1, 6
    		k = V2T(n,1)
    		l = V2T(n,2)
    		
    		C_tensor(i,j,k,l) = C_voigt(m,n) / (factor(m) * factor(n))
    		
    		! minor symmetries
    		C_tensor(j,i,k,l) = C_tensor(i,j,k,l)
    		C_tensor(i,j,l,k) = C_tensor(i,j,k,l)
    		C_tensor(j,i,l,k) = C_tensor(i,j,k,l)
    	END DO
    END DO
    
    RETURN
END SUBROUTINE VOIGT6x6_2_TENSOR3x3x3x3

    
!
! -------------------------------------------------------------------------------------------
! Tensor to Voigt transformation subroutine
! -------------------------------------------------------------------------------------------


SUBROUTINE TENSOR3x3x3x3_2_VOIGT6x6(C_tensor, C_voigt)

! Define input/output variables
! ----------------------------- 
    INTEGER :: i, j, k, l, m, n
    INTEGER :: V2T
    
    DOUBLE PRECISION :: C_voigt, C_tensor
    DOUBLE PRECISION :: factor

    DIMENSION C_voigt(6,6), C_tensor(3,3,3,3)
    DIMENSION factor(6), V2T(6,2)

! Float numbers / parameters
! --------------------------
    DOUBLE PRECISION :: Zero_R = 0.0D0, One_R = 1.0D0, Two_R = 2.0D0, Three_R = 3.0D0, Four_R = 4.0D0, Five_R = 5.0D0
    DOUBLE PRECISION :: Six_R = 6.0D0, Seven_R = 7.0D0, Eight_R = 8.0D0, Nine_R = 9.0D0, Ten_R = 10.0D0
    
! Assign values
! -------------
    ! For converting Engineering shear strain to true strain tensor
    factor(1) = One_R; factor(4) = Two_R
    factor(2) = One_R; factor(5) = Two_R
    factor(3) = One_R; factor(6) = Two_R
    
    ! Indices in Abaqus voigt notation
    V2T(1,1) = 1; V2T(1,2) = 1
    V2T(2,1) = 2; V2T(2,2) = 2
    V2T(3,1) = 3; V2T(3,2) = 3
    V2T(4,1) = 1; V2T(4,2) = 2
    V2T(5,1) = 1; V2T(5,2) = 3
    V2T(6,1) = 2; V2T(6,2) = 3
    
    ! Formulate the Stiffness tensor in Voigt notation
    DO m = 1, 6
    	i = V2T(m,1)
    	j = V2T(m,2)
    	DO n = 1, 6
    		k = V2T(n,1)
    		l = V2T(n,2)
    		
    		C_voigt(m,n) = C_tensor(i,j,k,l) * (factor(m) * factor(n))
    	END DO
    END DO
    
    RETURN
END SUBROUTINE TENSOR3x3x3x3_2_VOIGT6x6

