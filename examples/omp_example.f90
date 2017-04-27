program test_omp
!compile with "-fopenmp"
implicit none
integer, parameter:: n=1000000
integer a(n), b(n), c(n)
integer i
real start, fin

do i=1,n
        a(i)=i
        b(i)=2*i
end do

call cpu_time(start)
call system ('date')
        !#$OMP PARALLEL DO
        do i=1, n

        c(i)=myfun(a(i),b(i))
        end do
        !#$OMP END PARALLEL DO

call cpu_time(fin)
call system ('date')
        write(*,*) fin - start
write(*,*) sum(c)
contains

integer function myfun(dat1,dat2)
integer dat1,dat2,j, r
r = dat1
do j = 1,1000
        r = mod(r*dat2,165)
end do
myfun = r
end function

end program
