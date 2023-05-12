program main
use sistemas_lineales
    Implicit none
    interface
        function f(n,x)
        Implicit none
        integer, intent(in) :: n 
        real(8), intent(in) :: x(n) 
        real(8) :: f(n)
        end function f  

        function jf(n,x,h)
        Implicit none
        integer, intent(in) :: n 
        real(8), intent(in) :: x(n)
        real(8) :: jf(n,n)
        real(8), intent(in) :: h(n)
        end function jf 
    end interface 

    integer, parameter :: n = 3
    real(8) :: x0(n), B(n), h(n), x1(n), X(n), B0(n)
    integer :: i, kmax 
    real(8) :: A(n,n)
    real(8) :: ref, errorx, errory, eps, comprobar
    real(8), external :: norma
    
    eps = 1d-5
    kmax = 20
    x0 = (/0d0,0d0,0d0/)
    h = (/0.01d0,0.01d0,0.01d0/)

    B = -f(n,x0)
    !Resolución del método iterativo 
    do i = 1, kmax 
        A = jf(n,x0,h)
        
        call solve_gauss(3,A,B,X)!llamamos a la subrutina solve_gauss del módulo sistemas lineales para que resuelva el sistema
        
        x0 = x0 + X !x0 guarda el valor en el paso actual y lo utilizaremos para calcular el paso siguiente
        B = -f(n,x0)
    
        !Cálculo del error
        ref = norma(n,x0)
        errorx = norma(n,X) / ref !X guarda el vector de variación entre x0 y x1
        errory = norma(n,B)
        
        !Criterio de parada
        If(errorx < eps .and. errory < eps) then  
            write(*,*) "La solución al sistema de ecuaciones es :"
            Write(*,'(X,3(E16.10,X))') x0
        !comprobamos que la solución es correcta:
            B0 = f(n,x)
            If(B0(1)<1d-5 .and. B0(2)<1d-5 .and. B0(3)<1d-5) then
                Write(*,*) "La solución cumple la ecuación"
            end if  
                
        return
        end if  
    end do 
    write(*,*) "Cuidado: numero máximo de iteraciones alcanzado."
    stop




end program main
    function norma(n, x)
        !Calcula la norma de un vector
        implicit none
        integer, intent(in) :: n
        real(8), intent(in), dimension(n) :: x
        real(8) :: norma
        integer :: i
        norma = 0d0
        do i = 1, n
            norma = norma + x(i)**2
        end do
        norma = sqrt(norma)
    end function norma

    function f(n,x)
        Implicit none
        integer, intent(in) :: n 
        real(8), intent(in) :: x(n) 
        real(8) :: f(n)
        
        f(1) = (x(1)**2d0) + x(2) - 37d0
        f(2) = x(1) - (x(2)**2d0) - 5d0
        f(3) = x(1) + x(2) + x(3) - 3d0 
    end function f

    function jf(n,x,h)
        Implicit none
        integer, intent(in) :: n 
        real(8), intent(in) :: x(n)
        real(8) :: jf(n,n)
        real(8), intent(in) :: h(n)
    
        jf(1,1) = ((((x(1) + h(1))**2d0) + x(2) - 37d0) -  (((x(1) - h(1))**2d0) + x(2) - 37d0)) / (2d0*h(1))
        jf(1,2) = (((x(1)**2d0) + (x(2) + h(2)) - 37d0) - ((x(1)**2d0) + (x(2) - h(2)) - 37d0)) / (2d0*h(1))
        jf(1,3) = 0d0 

        jf(2,1) = (((x(1) + h(1)) - (x(2)**2d0) - 5d0) - ((x(1) - h(1)) - (x(2)**2d0) - 5d0)) / (2d0*h(2))
        jf(2,2) = ((x(1) - ((x(2) + h(2))**2d0) - 5d0) - (x(1) - ((x(2) - h(2))**2d0) - 5d0)) / (2d0*h(2))
        jf(2,3) = 0d0 

        jf(3,1) = (((x(1) + h(1)) + x(2) + x(3) - 3d0) - ((x(1) - h(1)) + x(2) + x(3) - 3d0)) / (2d0*h(3))
        jf(3,2) = ((x(1) + (x(2) + h(2)) + x(3) - 3d0) - (x(1) + (x(2) - h(2)) + x(3) - 3d0)) / (2d0*h(3))
        jf(3,3) = ((x(1) + x(2) + (x(3) + h(3)) - 3d0) - (x(1) + x(2) + (x(3) - h(3)) - 3d0)) / (2d0*h(3))

    end function jf

            


 