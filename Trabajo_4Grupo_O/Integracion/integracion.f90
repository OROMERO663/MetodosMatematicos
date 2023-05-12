module integracion
    implicit none
    contains
        subroutine trapecio(a,b,f,n,t)
        !calcula la integral entre a y b de f(x) con n subintervalos
            implicit none 
            integer :: n, i   
            real(8) :: a, b, x, h, t, s 
            real(8), external :: f 

            h = (b-a)/dble(n)
            !calculamos la suma de las ordenadas de 1 hasta n-1
            s = 0d0
            do i = 1, n-1
                x = a + dble(i)*h 
                s = s + f(x)
            end do 
            !fórmula del trapecio
            t = h * (0.5d0 * (f(a) + f(b)) + s)
            return 
        end subroutine trapecio
        
        subroutine trapecio_h(a,b,f,h,t)
            implicit none 
            integer :: n, i  
            real(8) :: a, b, x, h, t, s 
            real(8), external :: f 

            n = ceiling((b-a)/h)

            s = 0d0
            do i = 1, n-1
                x = a + dble(i)*h 
                s = s + f(x)
            end do 
            !fórmula del trapecio
            t = h * (0.5d0 * (f(a) + f(b)) + s)
            return 
        end subroutine trapecio_h  

        subroutine simpson_dup(a,b,f,eps,kmax,s)
            Implicit none
            integer :: n, kmax, k, i  
            real(8) :: a, b, t0, t1, eps, error, h, s, x, s0, sum  
            real(8), external :: f  
            !Calculamos el trapecio con 1 intervalo
            n = 1
            h = (b-a)/dble(n)
            t0 = 0.5d0 * h * (f(a) + f(b))
            !Calculamos el trapecio y simpson con dos subintervalos
            n = 2*n 
            h = 0.5*h 
            t1 = 0.5d0 * t0 + h * f(a + h)
            s = (4d0 * t1 - t0)/3d0
            !Guardamos t1 y s en t0 y s0 
            t0 = t1
            s0 = s 
            
            !Comenzamos las duplicaciones, n empieza siendo 4
            do k = 2, kmax
                n = 2**k
                h = (b-a)/dble(n) 
                !Calculamos el trapecio con n subintervalos
                !y calculamos la suma de las ordenadas impares
                sum = 0d0
                do i = 1, n/2 
                    x = a + dble(2*i-1)*h 
                    sum = sum + f(x)
                end do 
                !calculamos el nuevo trapecio
                t1 = 0.5d0 * t0 + h * sum 
                !calcula simpson con n subintervalos
                s = (4d0*t1 - t0)/3d0 
                !estimación del error
                error = abs((s-s0)/s0)
                !criterio de parada
                if (error < eps) then 
                return 
                end if 

                 write(*,'(I4,2X,2(E16.10,2X))') n, s, error 

                !Guardamos t1 y s en t0 y s0
                t0 = t1 
                s0 = s 
            end do 
            write(*,*) 'máximo número de duplicaciones'
            return 
        end subroutine simpson_dup

        function f(x)
        implicit none
        real(8) :: x, f  
            f = exp(x) * cos(x) 
        return 
        end function f
end module  

        !para el caso de que te den h en vez de n 
        !sacamos n = ceiling ((a-b)/h)
        !con ello sacamos el nuevo h y el nuevo n 
!ceiling(x) nos da el entero inmediatamente superior