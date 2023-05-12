! Metodos para la solucion de sistemas de ecuaciones lineales.
module sistemas_lineales

    implicit none
    contains

    subroutine solve_triang_inf(n,L,B,X)
        ! Resuelve el sistema de n ecuaciones triangular inferior definido por LX = B
        implicit none
        !Declaracion de argumentos de entrada:
        integer :: n  ! Tamaño del sistema de ecuaciones
        real(8), dimension(n,n) :: L  ! Matriz de coeficientes
        real(8), dimension(n) :: B  ! Vector de terminos independientes
        !Declaracion de argumentos de salida:
        real(8), dimension(n) :: X  ! Solucion del sistema de ecuaciones

        !Declaracion de variables auxiliares:
        integer :: k, h  ! Contadores
        real(8) :: s

        do k = 1, n  ! Resuelve las ecuaciones de la primera a la ultima
            ! Calculo de la suma
            s = 0d0
            do h = 1, k-1
                s = s + L(k,h)*X(h)
            end do
            ! Calculo de la solucion de la ecuacion
            X(k) = (B(k) - s)/L(k,k)
        end do
    end subroutine

    subroutine solve_triang_sup(n,U,B,X)
        ! Resuelve el sistema de n ecuaciones triangular superior definido por UX = B
        implicit none
        !Declaracion de argumentos de entrada:
        integer :: n  ! Tamaño del sistema de ecuaciones
        real(8), dimension(n,n) :: U  ! Matriz de coeficientes
        real(8), dimension(n) :: B  ! Vector de terminos independientes
        !Declaraci�n de argumentos de salida:
        real(8), dimension(n) :: X  ! Solucion del sistema de ecuaciones

        !Declaracion de variables auxiliares:
        integer :: k, h  ! Contadores
        real(8) :: s

        do k = n, 1,-1 ! Resuelve las ecuaciones de la ultima a la primera
            !Calculo de la suma
            s = 0d0
            do h = k+1, n
                s = s + U(k,h)*X(h)
            end do
            ! Resuelve la ecuacion
            X(k) = (B(k) - s)/U(k,k)
        end do
    end subroutine

    subroutine triangularizacion_gauss(n,A,B)
        ! Triangulariza el sistema de ecuaciones AX=B mediante operaciones fundamentales por filas
        ! utilizando el metodo de Gauss. La matriz de coeficientes triangularizada se guarda en la
        ! propia variable A y el vector de terminos independientes resultantes en B.
        implicit none
        ! Declaracion de argumentos de entrada y salida
        integer, intent(in) :: n  ! Tamaño del sistema de ecuaciones
        real(8), dimension(n,n), intent(inout) :: A  ! Matriz de coeficientes
        real(8), dimension(n), intent(inout) :: B  ! Matriz columna de terminos independientes.
        ! Declaracion de variables auxiliares
        integer :: k,i,j,r  ! Contadores
        real(8) :: c

        ! Avanza por la diagonal. En cada etapa hace ceros los elementos por debajo de la diagonal
        ! de la columna correspondiente.
        do k = 1, n-1

            call pivot(n,k,A,b) !aqui se llama a pivot
            do i = k+1,n ! Recorre las filas por debajo de k en orden descendente
                c = A(i,k)/A(k,k)  ! Multiplicador: Coeficiente que utiliza para modificar todos los elementos de la fila i
                do j=k,n  ! Recorre las columnas aplicando la operacion correspondiente
                    A(i,j) = A(i,j) - c*A(k,j)  ! Aplica la operacion fundamental
                end do
                B(i) = B(i) - c*B(k)  ! Aplica la misma operacion fundamental al vector de trminos ind.
            end do
            ! Para imprimir la matriz ampliada despues de cada etapa descomentar las siguientes 5 lineas.
            !write(*,*) "Matriz ampliada despues de etapa: ", k
            !do r = 1, n
            !    write(*, '(5(E13.6,XXX))') A(r,:), B(r)
            !end do
            !write(*,*) "---------------"
        end do
    end subroutine

    subroutine solve_gauss(n,A,B,X)
        !Resuelve el sistema de ecuaciones lineales AX=B mediante el metodo de Gauss
        implicit none
        ! Declaracion de los argumentos de entrada y salida
        integer :: n
        real(8), intent(inout), dimension(n,n) :: A  ! Matriz de coeficientes. Cuidado porque la matriz A se vera modificada en el proceso.
        real(8), intent(inout), dimension(n) :: B  ! Matriz columna de terminos independientes. Cuidado porque B se vera modificada en el proceso.
        real(8), intent(out), dimension(n) :: X  ! Solucion del sistema de ecuaciones.

        call triangularizacion_gauss(n, A, B)  ! Triangulariza la matriz ampliada
        call solve_triang_sup(n,A,B,X)  ! Resuelve el sistema triangular superior.
    end subroutine

    subroutine pivot(n,k,A,b)
        !selecciona el mayor elemento de cada columna para ponerlo como pivote a cada paso
        Implicit none
        integer, intent(in) :: n, k 
        real(8), intent(inout), dimension(n,n) :: A  
        real(8), intent(inout), dimension(n) :: b
        integer :: ii, jj, p 
        real(8) :: big, dummy 

        p = k !columna en la que se encuentra k
        big = abs(A(k,k)) !pivote actual que compararemos con los demás elementos

        do ii = k+1, n 
            dummy = abs(A(ii,k))
            if(dummy > big) then !Guardamos el elemento mayor y su fila asociadaen las nuevas variables temporales
                big = dummy 
                p = ii 
            end if  
        end do 

        if(p /= k) then 
            do jj = k, n !intercambiamos la fila del elemento original con la del pivote más grande
                dummy = A(p,jj)
                A(p,jj) = A(k,jj)
                A(k,jj) = dummy 
            end do 
            dummy = b(p) !hacemos que las filas de B sigan a las filas que estamos moviendo en A
            b(p) = b(k)
            b(k) = dummy 
        end if 
    end subroutine  

    subroutine solve_LU(n,A,B,X,Y,L,U)
        !resuelve el sistema de ecs AX=B por descomposicion LU
        Implicit none 
        integer :: n 
        real(8), intent(inout), dimension(n,n) :: A, L, U
        real(8), intent(inout), dimension(n) :: B, Y
        real(8), intent(out), dimension(n) :: X 
        call factor_LU(n,A,L,U)
        call solve_triang_infLU(n,L,B,Y)
        call solve_triang_supLU(n,U,Y,X)
    end subroutine 

    subroutine solve_triang_infLU(n,L,B,Y)
        ! Resuelve el sistema de n ecuaciones triangular inferior definido por LX = B
        implicit none
        !Declaracion de argumentos de entrada:
        integer, intent(in) :: n  ! Tamaño del sistema de ecuaciones
        real(8), intent(in), dimension(n,n) :: L  ! Matriz de coeficientes
        real(8), intent(in), dimension(n) :: B  ! Vector de terminos independientes
        !Declaracion de argumentos de salida:
        real(8), intent(out), dimension(n) :: Y  ! Solucion del sistema de ecuaciones

        !Declaracion de variables auxiliares:
        integer :: k, h  ! Contadores
        real(8) :: s

        do k = 1, n  ! Resuelve las ecuaciones de la primera a la ultima
            ! Calculo de la suma
            s = 0d0
            do h = 1, k-1
                s = s + L(k,h)*Y(h)
            end do
            ! Calculo de la solucion de la ecuacion
            Y(k) = (B(k) - s)/L(k,k)
        end do
    end subroutine

    subroutine solve_triang_supLU(n,U,Y,X)
        ! Resuelve el sistema de n ecuaciones triangular superior definido por UX = B
        implicit none
        !Declaracion de argumentos de entrada:
        integer, intent(in) :: n  ! Tamaño del sistema de ecuaciones
        real(8), intent(in), dimension(n,n) :: U  ! Matriz de coeficientes
        real(8), intent(in), dimension(n) :: Y  ! Vector de terminos independientes
        !Declaraci�n de argumentos de salida:
        real(8), intent(out), dimension(n) :: X  ! Solucion del sistema de ecuaciones

        !Declaracion de variables auxiliares:
        integer :: k, h  ! Contadores
        real(8) :: s

        do k = n, 1,-1 ! Resuelve las ecuaciones de la ultima a la primera
            !Calculo de la suma
            s = 0d0
            do h = k+1, n
                s = s + U(k,h)*X(h)
            end do
            ! Resuelve la ecuacion
            X(k) = (Y(k) - s)/U(k,k)
        end do
    end subroutine

    subroutine factor_LU(n,A,L,U) 
    !factoriza la matriz A en las matrices L y U 
    Implicit none 
    integer, intent(in) :: n 
    real(8), intent(in), dimension(n,n) :: A  
    real(8), intent(out), dimension(n,n) :: L, U 

    !Variables auxiliares del método
    integer :: i, k, j, h      
    real(8) :: multcol, multrow

    do k = 1, n !con este bucle avanzamos por la diagonal (k,k) a cada paso
        do i = k, n !bucle para calcular la columna k de L(i,k)
            multcol=0d0
            do h = 1, k-1
                
                multcol = multcol + (L(i,h) * U(h,k))
            end do 
            L(i,k) = A(i,k) - multcol   
        end do 

        do j = k+1, n !bucle para calcular la fila k de U(k,j) 
            multrow=0d0
            do h = 1, k-1
                
                multrow = multrow + (L(k,h) * U(h,j))  
            end do 
            U(k,j) = (1 / L(k,k)) * (A(k,j) - multrow)  
        end do 
    end do 
    !Bucles para generar los unos y ceros en las componentes restantes de L y U 
    do k = 1, n 
        U(k,k) = 1d0
    end do 
    do k = 2, n 
        do j = 1, k-1
            U(k,j) = 0d0
        end do 
    end do  
    do k = 2, n 
        do i = 1, k-1
            L(i,k) = 0d0
        end do
    end do 

    end subroutine

    subroutine solve_jacobi(n,A,B,X0,nmax,eps,X)
        !Resuelve el sistema de ecuaciones lineales AX=B mediante el metodo de Jacobi.
        implicit none
        integer, intent(in) :: n  ! Tamaño del sistema
        real(8), dimension(n,n), intent(in) :: A  ! Matriz de coeficientes
        real(8), dimension(n), intent(in) :: B  ! Vector de términos independientes
        real(8), dimension(n), intent(in) :: X0  ! Aproximacion inicial
        integer, intent(in) :: nmax  ! Numero maximo de iteraciones
        real(8), intent(in) :: eps  ! Tolerancia
        real(8), intent(out), dimension(n) :: X  ! Solución

        integer :: i, j, k
        real(8), dimension(n) :: X1, X2
        real(8) :: s1, s2
        real(8) :: e, ref

        X1 = X0  ! Inicializo la iteración anterior.

        do k = 1, nmax
            ! Aplicar ley de recurrencia para calcular la iteración siguiente
            do i = 1, n  ! Calcular las componentes una a una
                s1=0d0
                do j = 1, i-1
                    s1 = s1 + A(i,j)*X1(j)
                end do
                s2 =0d0
                do j = i+1, n
                    s2 = s2 + A(i,j)*X1(j)
                end do
                X2(i) = (B(i)-s1-s2)/A(i,i)
            end do

            ! Estimar el error
            e = norma2(X2-X1)
            ref = norma2(X2)
            if (ref>1d0) then
                e = e / ref
            end if

            ! Descomentar para imprimir iteraciones por pantalla.
        
            write(*,'(I2,X,4(E16.10,X))') k, X2, e

            ! Aplicar criterio de parada
            if (e<eps) then
                X = X2
                return
            end if
            X1 = X2  ! Actualizo el valor de la iteración anterior
        end do
        write(*,*) "Cuidado: numero máximo de iteraciones alcanzado."
        stop
    end subroutine

    subroutine solve_seidel(n,A,B,X0,nmax,eps,X)
        !Resuelve el sistema de ecuaciones lineales AX=B mediante el metodo de Jacobi.
        implicit none
        integer, intent(in) :: n  ! Tamaño del sistema
        real(8), dimension(n,n), intent(in) :: A  ! Matriz de coeficientes
        real(8), dimension(n), intent(in) :: B  ! Vector de términos independientes
        real(8), dimension(n), intent(in) :: X0  ! Aproximacion inicial
        integer, intent(in) :: nmax  ! Numero maximo de iteraciones
        real(8), intent(in) :: eps  ! Tolerancia
        real(8), intent(out), dimension(n) :: X  ! Solución

        integer :: i, j, k
        real(8), dimension(n) :: X1, X2
        real(8) :: s1, s2
        real(8) :: e, ref, eabs

        X1 = X0  ! Inicializo la iteración anterior
        do k = 1, nmax
            ! Aplicar ley de recurrencia para calcular la iteración siguiente
            do i = 1, n  ! Calcular las componentes una a una
                s1=0d0
                do j = 1, i-1
                    s1 = s1 + A(i,j)*X2(j)
                end do
                s2 =0d0
                do j = i+1, n
                    s2 = s2 + A(i,j)*X1(j)
                end do
                X2(i) = (B(i)-s1-s2)/A(i,i)
            end do

            ! Estimar el error
            eabs = norma2(X2-X1)
            e = norma2(X2-X1)
            ref = norma2(X2)
            if (ref>1d0) then
                e = e / ref
            end if

            ! Descomentar para imprimir iteraciones por pantalla.
            !ojo con el formato, si el vector X2 tiene mas de 4 componentes hay que cambiar el formato de escritura

            write(*,'(I2,X,6(E16.10,X))') k, X2, eabs, e

            ! Aplicar criterio de parada
            if (e<eps) then
                X = X2
                return
            end if
            X1 = X2  ! Actualizo el valor de la iteración anterior
        end do
        write(*,*) "Cuidado: numero máximo de iteraciones alcanzado."
        stop
    end subroutine

    subroutine inversa(n,AA,AI)
        !hay que hacer AX = I e ir resolviendo n sistemas de 
        !ecuaciones con gauss con pivote
        Implicit none
        !argumentos de entrada
        integer, intent(in) :: n  
        real(8), intent(inout), dimension(n,n) :: AA
        !argumento de salida
        real(8), intent(out), dimension(n,n) :: AI 
        !variables auxiliares
        integer :: i, j 
        real(8), dimension(n) :: B, X 
        real(8), dimension(n,n) :: A     
        B = 0d0 
        do i = 1, n !bucle para recorrer las columnas de la matriz inversa
            do j = 1, n  
                if(i .eq. j) then !con esto hacemos las componentes de B unos o ceros según toque
                    B(j) = 1d0
                else
                    B(j) = 0d0
                end if  
            end do
            A = AA !hay que reiniciar la matriz A a la original AA porque el metodo de gauss la va sobreescribiendo de un paso a otro
            call solve_gauss(n,A,B,X)

            do j = 1, n!con este bucle vamos guardando las componentes de la matriz inversa
                AI(j,i) = X(j)
            end do 
        end do 
        do i = 1, n 
            write(*, '(5(E13.6,XXX))') AI(i,:)
        end do 

    end subroutine

    function norma2(X)
        implicit none
        real(8), intent(in), dimension(:) :: X
        real(8) :: norma2
        norma2 = sqrt(dot_product(X,X))
    end function norma2

    function determinante(n,X)
        Implicit none
        real(8), intent(in), dimension(:,:) :: X 
        integer, intent(in) :: n 
        real(8) :: determinante
        integer :: i 
        determinante = 1d0
        do i = 1, n 
            determinante = determinante * X(i,i)
        end do 
    end function determinante 

    




end module

