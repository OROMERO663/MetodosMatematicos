program main
!Resuelve ejemplos con el método LU
use sistemas_lineales
Implicit none
real(8), allocatable, dimension(:,:) :: A, L, U 
real(8), allocatable, dimension(:) :: B, X, Y

write(*,*) "MÉTODO DE DESCOMPOSICÓN LU"
! Ejemplo 1:
write(*,*) "EJEMPLO 1"
write(*,*) "---------"
allocate(A(3,3),L(3,3),U(3,3),B(3),X(3),Y(3))  !meter n aqui
A = reshape((/3d0,-1d0,2d0,6d0,-1d0,4d0,9d0,-4d0,5d0/), (/3,3/))
B = (/0d0,2d0,1d0/)

call solve_LU(3,A,B,X,Y,L,U)

write(*,*) "La solución es:"
    
write(*, '(5(E13.6,XXX))') X

!Ejemplo 2:
write(*,*) "EJEMPLO 2"
write(*,*) "---------"
deallocate(A,B,X,L,U,Y)
allocate(A(4,4),L(4,4),U(4,4),B(4),X(4),Y(4))
A = reshape((/3d0,2d0,4d0,3d0,2d0,5d0,-1d0,1d0,&
0d0,-1d0,-8d0,-1d0,1d0,1d0,2d0,6d0/), (/4,4/))
B = (/1d0,2d0,-1d0,0d0/)

call solve_LU(4,A,B,X,Y,L,U)

write(*,*) "La solución es:"
    
write(*, '(5(E13.6,XXX))') X
end program main