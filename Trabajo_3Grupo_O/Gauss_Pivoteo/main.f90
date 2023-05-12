program main
    !Ejemplos para el método Gauss con pivoteo
    use sistemas_lineales
    implicit none
    real(8), allocatable, dimension(:,:) :: A
    real(8), allocatable, dimension(:) :: B, X

    write(*,*) "MÉTODO DE GAUSS CON PIVOTEO"
    ! Ejemplo 1:
    write(*,*) "EJEMPLO 1"
    write(*,*) "---------"
    allocate(A(3,3),B(3),X(3))  !meter n aqui
    A = reshape((/1d0,2d0,0d0,2d0,4d0,1d0,-1d0,0d0,-1d0/), (/3,3/))
    B = (/0d0,2d0,0d0/)

    call solve_gauss(3,A,B,X)  !meter n aqui

    write(*,*) "La solución es:"
    
    write(*, '(5(E13.6,XXX))') X

    !Ejemplo 2:
    write(*,*) "EJEMPLO 2"
    write(*,*) "---------"
    deallocate(A,B,X)
    allocate(A(4,4),B(4),X(4))
    A = reshape((/0d0,0d0,1d0,1d0,1d0,1d0,2d0,1d0,&
    1d0,1d0,-1d0,2d0,2d0,-1d0,3d0,0d0/), (/4,4/))
    B = (/2d0,-1d0,5d0,-2d0/)

    call solve_gauss(4,A,B,X)  !meter n aqui

    write(*,*) "La solución es:"
    
    write(*, '(5(E13.6,XXX))') X
    
end program main