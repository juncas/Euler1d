    subroutine tridiagsolve(n,a,b,c,d,x)
    implicit none
    integer, intent(in):: n                  ! dimension
    double precision, dimension(n):: a,b,c,d            !a(i)*x(i-1)+b(i)*x(i)+c(i)*x(i+1)=d(i)
    double precision, dimension(n),intent(out):: x 
    integer:: i       
    c(1)=c(1)/b(1)
    d(1)=d(1)/b(1)
    do i=2,n
        c(i)=c(i)/(b(i)-c(i-1)*a(i))
        d(i)=(d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i))
    end do
    x(n)=d(n)
    do i=n-1,1,-1
        x(i)=d(i)-c(i)*x(i+1)
    end do
    end