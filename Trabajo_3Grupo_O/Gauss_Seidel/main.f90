program main
!Ejemplos para el método de Gauss-Seidel
use sistemas_lineales
Implicit none
real(8), allocatable, dimension(:,:) :: A
real(8), allocatable, dimension(:) :: B, X
character(1), parameter :: L1 = "k"
character(3), parameter :: L2 = "X_1"
character(3), parameter :: L3 = "X_2"
character(3), parameter :: L4 = "X_3"
character(3), parameter :: L5 = "X_4"
character(9), parameter :: L6 = "error_abs"
character(9), parameter :: L7 = "error_rel"
write(*,*) "MÉTODO DE GAUSS-SEIDEL"

! Ejemplo 1:
!La aproximación inicial es X0 = (10,9,10)
!El número de iteraciones máxima es 50
!El error relativo suficiente para el criterio de parada es 1E-8
write(*,*) "EJEMPLO 1"
write(*,*) "---------"
allocate(A(3,3),B(3),X(3))  !meter n aqui
A = reshape((/5d0,1d0,2d0,2d0,6d0,1d0,-1d0,-3d0,4d0/), (/3,3/))
B = (/6d0,4d0,7d0/)
write(*,'(A1,X,3(7X,A3,8X),2(4X,A9,5X))') L1, L2, L3, L4, L6, L7
call solve_seidel(3,A,B,(/10d0,9d0,10d0/),50,1d-8,X)

write(*,*) "La solución es:"
write(*, '(X,3(E16.10,X))') X

!Ejemplo 2:
!La aproximación inicial es X0 = (10,9,10,10)
!El número de iteraciones máxima es 50
!El error relativo suficiente para el criterio de parada es 1E-8
write(*,*) "MÉTODO DE GAUSS-SEIDEL"
write(*,*) "EJEMPLO 2"
write(*,*) "---------"
deallocate(A,B,X)
allocate(A(4,4),B(4),X(4))  !meter n aqui
A = reshape((/5d0,1d0,2d0,1d0,2d0,6d0,1d0,-1d0,&
-1d0,-3d0,8d0,2d0,1d0,-3d0,-1d0,7d0/), (/4,4/))
B = (/6d0,4d0,7d0,2d0/)
write(*,'(A1,X,4(7X,A3,8X),2(4X,A9,5X))') L1, L2, L3, L4, L5, L6, L7
call solve_seidel(4,A,B,(/10d0,9d0,10d0,10d0/),50,1d-8,X)

write(*,*) "La solución es:"
write(*, '(X,4(E16.10,X))') X

end program main