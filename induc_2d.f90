! ====================================================================================================================

!	This program uses subroutine 'cyclic_trid' from the file 'trid_solve.f90'. 
!	Subroutine 'cyclic_trid' takes arguments as follows:
!	-----------------------------------------------------------------------------------------------------------------
!	cyclic_trid(dim, diag, sub, super, vec, X)
!	-----------------------------------------------------------------------------------------------------------------
!   PURPOSE:    solves a system of linear equations LX = D, where L is a cyclic tridigonal matrix
!               (i.e., tridiagonal + upper right and lower left corner elements are non-zero)  
!	-----------------------------------------------------------------------------------------------------------------
!	ARGUMENTS:  (in)	dim = dimensions of the system / number of unknowns
!               (in)	diag = vector containing diagonal elements of the matrix L (dimension dim)
!               (in)	sub = vector containing subdiagonal elements of the matrix L (dimension dim)
!                       		(the element sub(1) will be the upper right corner element of L)
!               (in)	super = vector containing superdiagonal elements of the matrix L (dimension dim)
!                       		(the element super(dim) will be the lower left corner element of L)
!               (in)	vec = the vector D in the above equation (dimension dim)
!               (out)	X = the vector containing the solutions (dimension dim)
!	-----------------------------------------------------------------------------------------------------------------

! =====================================================================================================================

! 	To run the program, write the following commands in terminal:
!		$ cd <file directory>
!		$ gfortran -o induc_2d induc_2d.f90 trid_solve.f90
!		$ ./induc_2d
!	The data of the field T(t, x, y) is stored in the file as specified near the bottom of the code.

! =====================================================================================================================

!	Before running the program, check the names of the output data files and see if they clash with any other 
!	filenames below. 
! 
!	-----------------------------------------------------------------------------------------------------------------
!	filename_1				filename_2/3/4				Vx(x, y)			Vy(x, y)			Plot
!	-----------------------------------------------------------------------------------------------------------------
!	induc_2d_dat_1			induc_2d_Xvals_1/			5					1					induc_2d_plot_1.gif
!							induc_2d_Yvals-1/
!							induc_2d_eta_t_1
!	-----------------------------------------------------------------------------------------------------------------
!	induc_2d_dat_2			induc_2d_Xvals_2/			0.3x				0.1y				induc_2d_plot_2.gif
!							induc_2d_Yvals-2/
!							induc_2d_eta_t_2
!	-----------------------------------------------------------------------------------------------------------------

! =====================================================================================================================

! 	The outputs of this program can be plotted using the program 'induc_2d_plot.py'.
!	
!	Change the filenames from which data is being read, as indicated in the table above, and enter an appropriate
!	output .gif filename.

! =====================================================================================================================

PROGRAM induc_2d

!	Define the parameters, e.g., dimension, time steps etc.

	IMPLICIT none
	INTEGER :: i, j, n, ierror, unit_no_1, unit_no_2, unit_no_3, unit_no_4, stepsize
	REAL, PARAMETER :: Lx = 5, Ly = 5, &
					   dx = 0.2, dy = 0.2
	REAL, PARAMETER :: time = 10, dt = 0.02
	INTEGER, PARAMETER :: Nx = int(2*Lx/dx) + 1, Ny = int(2*Ly/dy) + 1, steps = int(time/dt) + 1
	REAL, DIMENSION(steps, Nx, Ny) :: T
	REAL, DIMENSION(Nx, Ny) :: Vx, Vy, T_mid
	REAL, DIMENSION(Nx) :: X, vec1, diag1, sub1, super1, sol1
	REAL, DIMENSION(Ny) :: Y, vec2, diag2, sub2, super2, sol2
	REAL, PARAMETER :: eta = 0.3							! Diffusivity
	CHARACTER(len = 20) :: filename_1, filename_2, filename_3, filename_4
	CHARACTER(len = 100) :: frmt_1, frmt_2, frmt_3, frmt_4
	
	X = [(-Lx + (i - 1)*dx, i = 1, Nx)]
	Y = [(-Ly + (j - 1)*dy, j = 1, Ny)]
	
!	The output file names are listed below, change them as required.
	
	filename_1 = "induc_2d_dat_1"
	filename_2 = "induc_2d_Xvals_1"
	filename_3 = "induc_2d_Yvals_1"
	filename_4 = "induc_2d_eta_t_1"	
	
!	The integer 'stepsize' determines after how many time steps the value of T(n, i, j) will be stored.
!	stepsize*dt determines the time stepsize in the output plot.
!	Example:	stepsize = 25, dt = 0.02, time = 10
!				=> time stepsize in the plot = 25*0.02 = 0.5,
!				   number of frames = 10/0.5 + 1 = 21.

	stepsize = 25

!	Define the velocity fields Vx(x, y), Vy(x, y) and the initial condition T(t = 0, x, y) 
	
	DO j = 1, Ny
		DO i = 1, Nx
			Vx(i, j) = 5									! x-component of velocity field
			Vy(i, j) = 1									! y-component of velocity field
			T(1, i, j) = 10 * EXP(-(X(i)**2 + Y(j)**2))		! Initial condition
		END DO
	END DO

!	Coefficients in the implicit scheme
	
	DO i = 1, Nx
		diag1(i) = 1 + eta * (dt/dx**2)
		sub1(i) = -eta * (dt/(2*dx**2))
		super1(i) = -eta * (dt/(2*dx**2))
	END DO
	
	DO j = 1, Ny
		diag2(j) = 1 + eta * (dt/dy**2)
		sub2(j) = -eta * (dt/(2*dy**2))
		super2(j) = -eta * (dt/(2*dy**2))
	END DO

!	Solution using explicit + alternating direction implicit scheme
	
	time_loop: DO n = 1, steps - 1
		substep1: DO j = 1, Ny
			IF (j == 1) THEN
				vec1(1) = T(n, 1, 1) &
						  - (dt/(4*dx)) * (Vx(2, 1)*T(n, 2, 1) - Vx(Nx, 1)*T(n, Nx, 1)) &
						  - (dt/(4*dy)) * (Vy(1, 2)*T(n, 1, 2) - Vy(1, Ny)*T(n, 1, Ny)) &
						  + eta * (dt/(2*dy**2)) * (T(n, 1, 2) + T(n, 1, Ny) - 2*T(n, 1, 1))
				DO i = 2, Nx - 1
					vec1(i) = T(n, i, 1) &
						  	  - (dt/(4*dx)) * (Vx(i + 1, 1)*T(n, i + 1, 1) - Vx(i - 1, 1)*T(n, i - 1, 1)) &
						  	  - (dt/(4*dy)) * (Vy(i, 2)*T(n, i, 2) - Vy(i, Ny)*T(n, i, Ny)) &
						  	  + eta * (dt/(2*dy**2)) * (T(n, i, 2) + T(n, i, Ny) - 2*T(n, i, 1))
				END DO
				vec1(Nx) = T(n, Nx, 1) &
						  - (dt/(4*dx)) * (Vx(1, 1)*T(n, 1, 1) - Vx(Nx - 1, 1)*T(n, Nx - 1, 1)) &
						  - (dt/(4*dy)) * (Vy(Nx, 2)*T(n, Nx, 2) - Vy(Nx, Ny)*T(n, Nx, Ny)) &
						  + eta * (dt/(2*dy**2)) * (T(n, Nx, 2) + T(n, Nx, Ny) - 2*T(n, Nx, 1))
				CALL cyclic_trid(Nx, diag1, sub1, super1, vec1, sol1)
				DO i = 1, Nx
					T_mid(i, 1) = sol1(i)
				END DO
			ELSE IF (j == Ny) THEN
				vec1(1) = T(n, 1, Ny) &
						  - (dt/(4*dx)) * (Vx(2, Ny)*T(n, 2, Ny) - Vx(Nx, Ny)*T(n, Nx, Ny)) &
						  - (dt/(4*dy)) * (Vy(1, 1)*T(n, 1, 1) - Vy(1, Ny - 1)*T(n, 1, Ny - 1)) &
						  + eta * (dt/(2*dy**2)) * (T(n, 1, 1) + T(n, 1, Ny - 1) - 2*T(n, 1, Ny))
				DO i = 2, Nx - 1
					vec1(i) = T(n, i, Ny) &
						  	  - (dt/(4*dx)) * (Vx(i + 1, Ny)*T(n, i + 1, Ny) - Vx(i - 1, Ny)*T(n, i - 1, Ny)) &
						  	  - (dt/(4*dy)) * (Vy(i, 1)*T(n, i, 1) - Vy(i, Ny - 1)*T(n, i, Ny - 1)) &
						  	  + eta * (dt/(2*dy**2)) * (T(n, i, 1) + T(n, i, Ny - 1) - 2*T(n, i, Ny))
				END DO
				vec1(Nx) = T(n, Nx, Ny) &
						  - (dt/(4*dx)) * (Vx(1, Ny)*T(n, 1, Ny) - Vx(Nx - 1, Ny)*T(n, Nx - 1, Ny)) &
						  - (dt/(4*dy)) * (Vy(Nx, 1)*T(n, Nx, 1) - Vy(Nx, Ny - 1)*T(n, Nx, Ny - 1)) &
						  + eta * (dt/(2*dy**2)) * (T(n, Nx, 1) + T(n, Nx, Ny - 1) - 2*T(n, Nx, Ny))
				CALL cyclic_trid(Nx, diag1, sub1, super1, vec1, sol1)
				DO i = 1, Nx
					T_mid(i, Ny) = sol1(i)
				END DO
			ELSE
				vec1(1) = T(n, 1, j) &
						  - (dt/(4*dx)) * (Vx(2, j)*T(n, 2, j) - Vx(Nx, j)*T(n, Nx, j)) &
						  - (dt/(4*dy)) * (Vy(1, j + 1)*T(n, 1, j + 1) - Vy(1, j - 1)*T(n, 1, j - 1)) &
						  + eta * (dt/(2*dy**2)) * (T(n, 1, j + 1) + T(n, 1, j - 1) - 2*T(n, 1, j))
				DO i = 2, Nx - 1
					vec1(i) = T(n, i, j) &
						  	  - (dt/(4*dx)) * (Vx(i + 1, j)*T(n, i + 1, j) - Vx(i - 1, j)*T(n, i - 1, j)) &
						  	  - (dt/(4*dy)) * (Vy(i, j + 1)*T(n, i, j + 1) - Vy(i, j - 1)*T(n, i, j - 1)) &
						  	  + eta * (dt/(2*dy**2)) * (T(n, i, j + 1) + T(n, i, j - 1) - 2*T(n, i, j))
				END DO
				vec1(Nx) = T(n, Nx, j) &
						  - (dt/(4*dx)) * (Vx(1, j)*T(n, 1, j) - Vx(Nx - 1, j)*T(n, Nx - 1, j)) &
						  - (dt/(4*dy)) * (Vy(Nx, j + 1)*T(n, Nx, j + 1) - Vy(Nx, j - 1)*T(n, Nx, j - 1)) &
						  + eta * (dt/(2*dy**2)) * (T(n, Nx, j + 1) + T(n, Nx, j - 1) - 2*T(n, Nx, j))
				CALL cyclic_trid(Nx, diag1, sub1, super1, vec1, sol1)
				DO i = 1, Nx
					T_mid(i, j) = sol1(i)
				END DO
			END IF
		END DO substep1
		
		substep2: DO i = 1, Nx
			IF (i == 1) THEN
				vec2(1) = T_mid(1, 1) &
						- (dt/(4*dx)) * (Vx(2, 1)*T_mid(2, 1) - Vx(Nx, 1)*T_mid(Nx, 1)) &
						- (dt/(4*dy)) * (Vy(1, 2)*T_mid(1, 2) - Vy(1, Ny)*T_mid(1, Ny)) &
						+ eta * (dt/(2*dx**2)) * (T_mid(2, 1) + T_mid(Nx, 1) - 2*T_mid(1, 1))
				DO j = 2, Ny - 1
					vec2(j) = T_mid(1, j) &
							  - (dt/(4*dx)) * (Vx(2, j)*T_mid(2, j) - Vx(Nx, j)*T_mid(Nx, j)) &
							  - (dt/(4*dy)) * (Vy(1, j + 1)*T_mid(1, j + 1) - Vy(1, j - 1)*T_mid(1, j - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(2, j) + T_mid(Nx, j) - 2*T_mid(1, j))
				END DO
				vec2(Ny) = T_mid(1, Ny) &
							  - (dt/(4*dx)) * (Vx(2, Ny)*T_mid(2, Ny) - Vx(Nx, Ny)*T_mid(Nx, Ny)) &
							  - (dt/(4*dy)) * (Vy(1, 1)*T_mid(1, 1) - Vy(1, Ny - 1)*T_mid(1, Ny - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(2, Ny) + T_mid(Nx - 1, Ny) - 2*T_mid(1, Ny))
				CALL cyclic_trid(Ny, diag2, sub2, super2, vec2, sol2)
				DO j = 1, Ny
					T(n + 1, 1, j) = sol2(j)
				END DO
			ELSE IF (i == Nx) THEN
				vec2(1) = T_mid(Nx, 1) &
						- (dt/(4*dx)) * (Vx(1, 1)*T_mid(1, 1) - Vx(Nx - 1, 1)*T_mid(Nx - 1, 1)) &
						- (dt/(4*dy)) * (Vy(Nx, 2)*T_mid(Nx, 2) - Vy(Nx, Ny)*T_mid(Nx, Ny)) &
						+ eta * (dt/(2*dx**2)) * (T_mid(1, 1) + T_mid(Nx - 1, 1) - 2*T_mid(Nx, 1))
				DO j = 2, Ny - 1
					vec2(j) = T_mid(Nx, j) &
							  - (dt/(4*dx)) * (Vx(1, j)*T_mid(1, j) - Vx(Nx - 1, j)*T_mid(Nx - 1, j)) &
							  - (dt/(4*dy)) * (Vy(Nx, j + 1)*T_mid(Nx, j + 1) - Vy(Nx, j - 1)*T_mid(Nx, j - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(1, j) + T_mid(Nx - 1, j) - 2*T_mid(Nx, j))
				END DO
				vec2(Ny) = T_mid(i, Ny) &
							  - (dt/(4*dx)) * (Vx(i + 1, Ny)*T_mid(i + 1, Ny) - Vx(i - 1, Ny)*T_mid(i - 1, Ny)) &
							  - (dt/(4*dy)) * (Vy(i, 1)*T_mid(i, 1) - Vy(i, Ny - 1)*T_mid(i, Ny - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(i + 1, Ny) + T_mid(i - 1, Ny) - 2*T_mid(i, Ny))
				CALL cyclic_trid(Ny, diag2, sub2, super2, vec2, sol2)
				DO j = 1, Ny
					T(n + 1, i, j) = sol2(j)
				END DO
			ELSE
				vec2(1) = T_mid(i, 1) &
						- (dt/(4*dx)) * (Vx(i + 1, 1)*T_mid(i + 1, 1) - Vx(i - 1, 1)*T_mid(i - 1, 1)) &
						- (dt/(4*dy)) * (Vy(i, 2)*T_mid(i, 2) - Vy(i, Ny)*T_mid(i, Ny)) &
						+ eta * (dt/(2*dx**2)) * (T_mid(i + 1, 1) + T_mid(i - 1, 1) - 2*T_mid(i, 1))
				DO j = 2, Ny - 1
					vec2(j) = T_mid(i, j) &
							  - (dt/(4*dx)) * (Vx(i + 1, j)*T_mid(i + 1, j) - Vx(i - 1, j)*T_mid(i - 1, j)) &
							  - (dt/(4*dy)) * (Vy(i, j + 1)*T_mid(i, j + 1) - Vy(i, j - 1)*T_mid(i, j - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(i + 1, j) + T_mid(i - 1, j) - 2*T_mid(i, j))
				END DO
				vec2(Ny) = T_mid(i, Ny) &
							  - (dt/(4*dx)) * (Vx(i + 1, Ny)*T_mid(i + 1, Ny) - Vx(i - 1, Ny)*T_mid(i - 1, Ny)) &
							  - (dt/(4*dy)) * (Vy(i, 1)*T_mid(i, 1) - Vy(i, Ny - 1)*T_mid(i, Ny - 1)) &
							  + eta * (dt/(2*dx**2)) * (T_mid(i + 1, Ny) + T_mid(i - 1, Ny) - 2*T_mid(i, Ny))
				CALL cyclic_trid(Ny, diag2, sub2, super2, vec2, sol2)
				DO j = 1, Ny
					T(n + 1, i, j) = sol2(j)
				END DO
			END IF
		END DO substep2
	END DO time_loop

!	Save the values of T(n, i, j) in a file
	
	unit_no_1 = 12
	
	WRITE (frmt_1, "(A, I0, A)") "(", Ny, "(1X, F5.2))"
	
	OPEN (UNIT = unit_no_1, FILE = filename_1, STATUS = 'UNKNOWN', ACTION = 'WRITE', IOSTAT = ierror)
	DO n = 1, steps, stepsize
		DO i = 1, Nx
			WRITE (unit_no_1, frmt_1) ((T(n, i, j)), j = 1, Ny)
		END DO
	END DO
	
	unit_no_2 = 13
	unit_no_3 = 14
	unit_no_4 = 15
	
	WRITE (frmt_2, "(A, I0, A)") "(", Nx, "(1X, F4.1))"
	WRITE (frmt_3, "(A, I0, A)") "(", Ny, "(1X, F4.1))"
	WRITE (frmt_4, "(A, I0, A)") "(", steps + 1, "(1X, F5.2))"
	
	OPEN(UNIT = unit_no_2, FILE = filename_2, STATUS = 'UNKNOWN', ACTION = 'WRITE', IOSTAT = ierror)
	WRITE (unit_no_2, frmt_2) ((X(i)), i = 1, Nx)
	
	OPEN(UNIT = unit_no_3, FILE = filename_3, STATUS = 'UNKNOWN', ACTION = 'WRITE', IOSTAT = ierror)
	WRITE (unit_no_3, frmt_3) ((Y(j)), j = 1, Ny)
	
	OPEN(UNIT = unit_no_4, FILE = filename_4, STATUS = 'UNKNOWN', ACTION = 'WRITE', IOSTAT = ierror)
	WRITE (unit_no_4, frmt_4) eta, (((n - 1)*dt), n = 1, steps, 25)

END PROGRAM induc_2d

! =====================================================================================================================
!	1.	Note that Lax scheme has not been implemented in the explicit part of the difference equation. This might
!		result in some instabilities for really bad choices of dx and dt, however, it seems satisfying the Courant
!		condition is sufficient for obtaining stable soultions.
























