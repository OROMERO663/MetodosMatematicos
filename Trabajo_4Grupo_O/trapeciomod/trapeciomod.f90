program main
use integracion
    Implicit none
    integer :: n, p, i
    real(8) :: a, b, t, t1, pi, h, e, eps, G 
    character(1), parameter :: L1 = "n"
    character(3), parameter :: L2 = "T_n"
    character(5), parameter :: L3 = "error"

    pi = acos(-1d0)
    G = -12.0703463164d0 !resultado analítico de la integral, lo vamos a usar para sacar el error
    p = 20 
    a = 0d0
    b = pi 
    eps = 1d-8

    write(*,*) "Método del Trapecio para tolerancia < 10**-8"
    write(*,'(2X,A1,11X,A3,11X,A5)') L1, L2, L3

    do i = 1, p !iteraciones para calcular el paso
        n = 2**i   
        h = (b-a)/dble(n)
        call trapecio_h(a,b,f,h,t)
        !Estimación del error
        e = abs((G-t)/G)
    
        write(*,'(I4,2X,2(E16.10,2X))') n, t, e 
    !criterio de parada
    If (e<eps) then 
        t = t1
        return 
        end if 
     
    end do
    stop 
    
        
end program main