program Ejercicio1
implicit none
integer :: n, i     
real(8), allocatable :: X(:), Y(:)

    open(15, file='test_sort.csv')
    read(15, '(19X,I2)') n !numero de componentes del vector
    write(*,*) n 

    allocate(X(n), Y(n))
    do i =1,n 
        read(15, '(3X, 2(XX,E23.16,X))') X(i), Y(i)
    end do 
    close(15)
    do i =1, n 
        write(*,*) X(i), Y(i)
    end do 
write(*,*) "-------------------------"
    call insertion_sort(n,X,Y)
    do i =1, n 
        write(*,*) X(i), Y(i)
    end do


end program Ejercicio1

subroutine insertion_sort(n,X,Y)
    Implicit none
    integer, intent(in) :: n 
    integer :: i, j 
    real(8) :: temp, tempy!variable donde guardaremos el número a evaluar/mover temporalmente
    real(8), dimension(n), intent(in out) :: X, Y
    !Sentencias ejecutables
     !inicialización de la variable

    do j = 2, n  
        temp = X(j)
        tempy = Y(j)
        i = (j - 1)
        do while ((i>0) .and. (X(i)>temp)) 
            X(i+1) = X(i)
            Y(i+1) = Y(i)
            i = (i - 1)
        end do
        X(i+1)=temp 
        Y(i+1)=tempy
    end do 

end subroutine insertion_sort 