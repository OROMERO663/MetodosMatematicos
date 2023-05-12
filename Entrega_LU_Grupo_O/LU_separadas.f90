program desc_lu
    implicit none 
    real(8), allocatable, dimension(:,:) :: A, L, U    
    integer :: i, k, j, h, r, n      
    real(8) :: multcol, multrow !los vamos a usar para los sumatorios 

    n = 3
    allocate(A(n,n), L(n,n), U(n,n))
    !Inicialización de la matriz
    A = reshape((/3d0,-1d0,2d0,6d0,-1d0,4d0,9d0,-4d0,5d0/),(/3,3/))

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
        
    !Escribimos las matrices por pantalla
    do r = 1, n
        write(*, '(5(E13.6,XXX))') A(r,:)
    end do 
    write(*,*) "------------------------"
    do r = 1, n
        write(*, '(5(E13.6,XXX))') L(r,:)
    end do 
    write(*,*) "------------------------"
    do r = 1, n
        write(*, '(5(E13.6,XXX))') U(r,:)
    end do 

    
end program desc_lu