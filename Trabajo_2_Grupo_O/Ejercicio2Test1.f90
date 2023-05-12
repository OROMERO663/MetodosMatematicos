program Ejecicio2
    implicit none
    integer :: n, i  
    real(8), allocatable :: X(:), Y(:)
    real(8) :: z
    logical :: sorted 
    real(8), external :: interp1d
    !Caso ejemplo 1
    open(16, file='test_interp_1.csv')
    read(16, '(19X, I2)') n !cantidad de puntos

    allocate(X(n), Y(n))

    do i = 1, n 
        read(16, '(3X, 2(XX,E23.16,X))') X(i), Y(i)
    end do 
    close(16)

    !comprobación de carga de puntos
    do i = 1, n 
        write(*,*) X(i), Y(i)
    end do  
        write(*,*) "-----------------------------------"
    !comprobamos que el vector X está ordenado
    sorted = .true.
    do i = 2, n
        if(X(i) <= X(i-1)) then
            sorted = .False.
            exit 
        end if
    end do


    z=2d0
    write(*,*) "para z = 2.0 el valor y(z) interpolado es: "
    write(*,*) interp1d(n,X,Y,z,sorted)

    z=5.5d0
    write(*,*) "para z = 5.5 el valor y(z) interpolado es: "
    write(*,*) interp1d(n,X,Y,z,sorted)

end program Ejecicio2
function interp1d(n,X,Y,z,sorted)
    implicit none
    integer, intent(in) :: n 
    real(8), intent(in), dimension(n) :: X, Y
    real(8), intent(in) :: z 
    logical, intent(in) :: sorted
    real(8) :: interp1d

    integer :: j, i  
    real(8) :: x1, x2, y1, y2 

if (sorted .eqv. .False.) then  
    call insertion_sort(n,X,Y)
end if  
 
    do i = 2, n 
        do j = 2, n  
            if ((X(j) >= z) .and. (X(j-1) < z)) then
                x1 = X(j-1)
                x2 = X(j)
                y1 = Y(j-1)
                y2 = Y(j)
            end if  
        end do
    end do 

    interp1d = y1 + ((y2 - y1)/(x2 - x1)) * (z -x1)
    
end function interp1d


