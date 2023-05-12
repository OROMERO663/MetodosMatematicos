program main
use integracion
    Implicit none
    real(8) :: a, b, s, pi, eps, error
    integer :: n, kmax 
    
    character(1), parameter :: L1 = "n"
    character(3), parameter :: L2 = "S_n"
    character(5), parameter :: L3 = "error"
    pi = acos(-1d0)
    write(*,*) "Método de Simpson para tolerancia < 10**-5"
    write(*,'(2X,A1,11X,A3,11X,A5)') L1, L2, L3
    !prueba para tolerancia menor de 10**-5
    a = 0d0
    b = pi 
    eps = 1d-5
    kmax = 20 

    call simpson_dup(a,b,f,eps,kmax,s)
    write(*,*) "-------------------------------------------"
    write(*,*) "Método de Simpson para tolerancia < 10**-8"
    write(*,'(2X,A1,11X,A3,11X,A5)') L1, L2, L3
    !prueba para tolerancia menor de 10**-8
    a = 0d0
    b = pi 
    eps = 1d-8
    kmax = 20 
    call simpson_dup(a,b,f,eps,kmax,s)

end program main 
