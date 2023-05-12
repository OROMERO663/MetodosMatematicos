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