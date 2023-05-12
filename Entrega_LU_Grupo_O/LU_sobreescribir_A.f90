program desc_lu
    implicit none 
    real(8), allocatable, dimension(:,:) :: A   
    integer :: i, k, j, h, r, n      
    real(8) :: multcol, multrow !los vamos a usar para los sumatorios 

    n = 3
    allocate(A(n,n)
    A = reshape((/3d0,-1d0,2d0,6d0,-1d0,4d0,9d0,-4d0,5d0/),(/3,3/))

    do k = 1, n !con este bucle avanzamos por la diagonal (k,k) a cada paso
        do i = k, n !bucle para calcular la columna k de L(i,k)
            multcol=0d0 
            do h = 1, k-1 
                multcol = multcol + (A(i,h) * A(h,k))
            end do 
            A(i,k) = A(i,k) - multcol   
        end do 

        do j = k+1, n !bucle para calcular la fila k de U(k,j) 
            multrow=0d0
            do h = 1, k-1               
                multrow = multrow + (A(k,h) * A(h,j))  
            end do 
            A(k,j) = (1 / A(k,k)) * (A(k,j) - multrow)  
        end do 
    end do 
!U y L están almacenados ahora en A
    do r = 1, n
        write(*, '(5(E13.6,XXX))') A(r,:)
    end do  
end program desc_lu