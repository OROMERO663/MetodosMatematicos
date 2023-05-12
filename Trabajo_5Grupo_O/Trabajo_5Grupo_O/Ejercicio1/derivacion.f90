program main
    Implicit none 
    integer :: n, i  
    real(8), parameter :: pi = acos(-1d0)
    real(8), parameter :: a = 0d0, b = 2*pi 
    real(8), allocatable, dimension(:) :: x , df, coseno  
    real(8) :: k, h, sum, MAE50, MAE100, MAE200, MAE500, MAE1000  !intervalos y sumatorio 
    real(8) :: f 

    n = 50
    allocate(x(n), df(n), coseno(n))

    k = (b-a)/dble(n)

    do i = 1, n 
        x(i) = a + dble(i+1) * k 
    end do 

    h = (x(3) - x(2))
    !diferencias adelantadas para el punto a
    df(1) = (-3d0*f(x(1)) + 4d0*f(x(1)+h) - f(x(1)+2d0*h)) / (2d0*h)
    !diferencias atrasadas para el punto b 
    df(n) = (3d0*f(x(n)) - 4d0*f(x(n)-h) + f(x(n)-2d0*h)) / (2d0*h)
    
    do i = 2, n-1 
        df(i) = ((f(x(i)+h) - f(x(i)-h))) / (2d0*h) 
    end do 
    do i = 1, n 
        coseno(i) = cos(x(i))
    end do 

    !Cálculo del error
    sum = 0d0
    do i = 1, n 
        sum = sum + abs(df(i) - coseno(i))
    end do 

    MAE50 = (1d0/dble(n+1)) * sum 

    open(15, file='derivación50.txt')
    write(15,'(A6,XX,E16.10,2X)') "MAE50:", MAE50
    write(15,'(A17,XX,I2)') "Numero de puntos:", n
    do i = 1, n 
        write(15,'(I2,2X,3(E16.10,2X))') i, x(i), df(i), coseno(i)   
    end do 
    close(15)

    deallocate(x, df, coseno)
    n = 100
    allocate(x(n), df(n), coseno(n))

    k = (b-a)/dble(n)

    do i = 1, n 
        x(i) = a + dble(i+1) * k 
    end do 

    h = (x(3) - x(2))
    !diferencias adelantadas para el punto a
    df(1) = (-3d0*f(x(1)) + 4d0*f(x(1)+h) - f(x(1)+2d0*h)) / (2d0*h)
    !diferencias atrasadas para el punto b 
    df(n) = (3d0*f(x(n)) - 4d0*f(x(n)-h) + f(x(n)-2d0*h)) / (2d0*h)
    
    do i = 2, n-1 
        df(i) = ((f(x(i)+h) - f(x(i)-h))) / (2d0*h) 
    end do 
    do i = 1, n 
        coseno(i) = cos(x(i))
    end do 

    !Cálculo del error
    sum = 0d0
    do i = 1, n 
        sum = sum + abs(df(i) - coseno(i))
    end do 

    MAE100 = (1d0/dble(n+1)) * sum 

    open(16, file='derivación100.txt')
    write(16,'(A6,XX,E16.10,2X)') "MAE100:", MAE100
    write(16,'(A17,XX,I3)') "Numero de puntos:", n
    do i = 1, n 
        write(16,'(I2,2X,3(E16.10,2X))') i, x(i), df(i), coseno(i)   
    end do 
    close(16)

    deallocate(x, df, coseno)
    n = 200
    allocate(x(n), df(n), coseno(n))

    k = (b-a)/dble(n)

    do i = 1, n 
        x(i) = a + dble(i+1) * k 
    end do 

    h = (x(3) - x(2))
    !diferencias adelantadas para el punto a
    df(1) = (-3d0*f(x(1)) + 4d0*f(x(1)+h) - f(x(1)+2d0*h)) / (2d0*h)
    !diferencias atrasadas para el punto b 
    df(n) = (3d0*f(x(n)) - 4d0*f(x(n)-h) + f(x(n)-2d0*h)) / (2d0*h)
    
    do i = 2, n-1 
        df(i) = ((f(x(i)+h) - f(x(i)-h))) / (2d0*h) 
    end do 
    do i = 1, n 
        coseno(i) = cos(x(i))
    end do 

    !Cálculo del error
    sum = 0d0
    do i = 1, n 
        sum = sum + abs(df(i) - coseno(i))
    end do 

    MAE200 = (1d0/dble(n+1)) * sum 

    open(17, file='derivación200.txt')
    write(17,'(A6,XX,E16.10,2X)') "MAE200:", MAE200
    write(17,'(A17,XX,I3)') "Numero de puntos:", n
    do i = 1, n 
        write(17,'(I3,2X,3(E16.10,2X))') i, x(i), df(i), coseno(i)   
    end do 
    close(17)

    deallocate(x, df, coseno)
    n = 500
    allocate(x(n), df(n), coseno(n))

    k = (b-a)/dble(n)

    do i = 1, n 
        x(i) = a + dble(i+1) * k 
    end do 

    h = (x(3) - x(2))
    !diferencias adelantadas para el punto a
    df(1) = (-3d0*f(x(1)) + 4d0*f(x(1)+h) - f(x(1)+2d0*h)) / (2d0*h)
    !diferencias atrasadas para el punto b 
    df(n) = (3d0*f(x(n)) - 4d0*f(x(n)-h) + f(x(n)-2d0*h)) / (2d0*h)
    
    do i = 2, n-1 
        df(i) = ((f(x(i)+h) - f(x(i)-h))) / (2d0*h) 
    end do 
    do i = 1, n 
        coseno(i) = cos(x(i))
    end do 

    !Cálculo del error
    sum = 0d0
    do i = 1, n 
        sum = sum + abs(df(i) - coseno(i))
    end do 

    MAE500 = (1d0/dble(n+1)) * sum 

    open(18, file='derivación500.txt')
    write(18,'(A6,XX,E16.10,2X)') "MAE500:", MAE500
    write(18,'(A17,XX,I3)') "Numero de puntos:", n
    do i = 1, n 
        write(18,'(I3,2X,3(E16.10,2X))') i, x(i), df(i), coseno(i)   
    end do 
    close(18)

    deallocate(x, df, coseno)
    n = 1000
    allocate(x(n), df(n), coseno(n))

    k = (b-a)/dble(n)

    do i = 1, n 
        x(i) = a + dble(i+1) * k 
    end do 

    h = (x(3) - x(2))
    !diferencias adelantadas para el punto a
    df(1) = (-3d0*f(x(1)) + 4d0*f(x(1)+h) - f(x(1)+2d0*h)) / (2d0*h)
    !diferencias atrasadas para el punto b 
    df(n) = (3d0*f(x(n)) - 4d0*f(x(n)-h) + f(x(n)-2d0*h)) / (2d0*h)
    
    do i = 2, n-1 
        df(i) = ((f(x(i)+h) - f(x(i)-h))) / (2d0*h) 
    end do 
    do i = 1, n 
        coseno(i) = cos(x(i))
    end do 

    !Cálculo del error
    sum = 0d0
    do i = 1, n 
        sum = sum + abs(df(i) - coseno(i))
    end do 

    MAE1000 = (1d0/dble(n+1)) * sum 

    open(19, file='derivación1000.txt')
    write(19,'(A6,XX,E16.10,2X)') "MAE1000:", MAE1000
    write(19,'(A17,XX,I4)') "Numero de puntos:", n
    do i = 1, n 
        write(19,'(I4,2X,3(E16.10,2X))') i, x(i), df(i), coseno(i)   
    end do 
    close(19)



end program main

function f(x)
    implicit none
    real(8) :: x, f  
        f = sin(x)
    return 
end function f