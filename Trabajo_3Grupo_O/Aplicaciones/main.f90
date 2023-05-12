program main
    use sistemas_lineales
    Implicit none
    integer, parameter :: n=4
    real(8), dimension(n,n) :: A, L, U 
    real(8), dimension(n,n) :: AA, AI
    real(8) :: det 

    A = reshape((/1d0,2d0,-1d0,1d0,4d0,3d0,-0.5d0,1d0,&
    1d0,1d0,2d0,-4d0,-3d0,-2d0,-2d0,1d0/),(/4,4/)) 
    AA = reshape((/1d0,2d0,-1d0,1d0,4d0,3d0,-0.5d0,1d0,&
    1d0,1d0,2d0,-4d0,-3d0,-2d0,-2d0,1d0/),(/4,4/)) 

    call factor_LU(n,A,L,U)
    det = determinante(n,L)

    write(*,*) "El determinante de la matriz A es: "
    write(*,*) det 
    write(*,*) "------------------"
    write(*,*) "La inversa de la matriz A es: "
    call inversa(n,AA,AI)

end program main 
