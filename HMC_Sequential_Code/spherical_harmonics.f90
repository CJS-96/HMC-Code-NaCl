module spherical_harmonics
  use consts, only : PR, PI
  implicit none
  private 
  public :: sp_hrmcs, wigner
contains
    !***************************************************************************
    !**          Computes the spherical harmonics Y(l,m,theta,phi)            **
    !***************************************************************************
    function sp_hrmcs(l,m,theta,phi)
        integer, intent(in) :: l,m
        real(PR), intent(in) :: theta,phi
        complex(PR) :: sp_hrmcs

        real(PR) :: z,fx,fy,fz
        integer :: i,mc

        fx=real(2*l+1,PR)/4._pr/pi
        if(l-m == 0)then
            fy=1._pr
        else
            fy=1._pr
            do i=1,l-m
                fy=fy*real(i,PR)
            enddo
        endif
        if(l+m == 0)then
            fz=1._pr
        else
            fz=1._pr
            do i=1,l+m
                fz=fz*real(i,PR)
            enddo
        endif
        z=cos(theta)

        if(m >= 0)then
            sp_hrmcs=sqrt(fx*fy/fz)*plgndr_s(l,m,z)*exp((0,1)*m*phi)
        else
            mc=abs(m)
            sp_hrmcs=sqrt(fx*fy/fz)*((-1)**mc*fz/fy*plgndr_s(l,mc,z))*exp((0,1)*m*phi)
        endif

    end function sp_hrmcs
    !****************************************************************************
    !** Computes the associated Legendre polynomial. Here l and m are integers **
    !** satisfying 0 <= m <= l, while x lies in the range -1 <= x <= 1.        **
    !**          (Ref: Numerical Recipies in Fortran 90, pg. 1122)             **
    !****************************************************************************
    function plgndr_s(l,m,x)
        integer, intent(in) :: l,m
        real(PR), intent(in) :: x
        real(PR) :: plgndr_s

        integer :: ll
        real(PR) :: pll,pmm,pmmp1,somx2

        if(.not.(m >= 0 .and. m <= l .and. abs(x) <= 1._pr))then
            write(*,*)'nrerror: an assertion failed with this tag: plgndr_s args'
            stop 'program terminated in plgndr_s'
        endif

        pmm=1._pr
        if(m > 0)then
            somx2=sqrt((1._pr-x)*(1._pr+x))
            pmm=product(arth(1._pr,2._pr,m))*somx2**m
            if(mod(m,2) == 1)pmm=-pmm
        endif
        if(l == m)then
            plgndr_s=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if( l == m+1)then
                plgndr_s=pmmp1
            else
                do ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                enddo
                plgndr_s=pll
            endif
        endif
    end function plgndr_s

    !****************************************************************************
    !Array function returning an arithmetic progression                        **
    !****************************************************************************
    pure function arth(first,increment,n)
        real(PR), intent(in) :: first,increment
        integer, intent(in) :: n
        real(PR), dimension(n) :: arth

        integer :: k, k2
        real(PR) :: temp

        integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8

        if(n > 0)arth(1)=first
        if(n <= NPAR_ARTH)then
            do k=2,n
                arth(k)=arth(k-1)+increment
            enddo
        else
            do k=2,NPAR2_ARTH
                arth(k)=arth(k-1)+increment
            enddo
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if(k >= n)exit
                k2=k+k
                arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                temp=temp+temp
                k=k2
            enddo
        endif
    end function arth

    !*****************************************************************************
    !** Function calculating wigner 3-j coefficients                            **
    !*****************************************************************************
    function wigner(l1,l2,l3)
        integer, intent(in) :: l1,l2,l3
        real(PR), dimension(-l1:l1,-l2:l2,-l3:l3) :: wigner

        integer :: m1,m2,m3

        wigner=0._pr
        do m1=-l1,l1
            do m2=-l2,l2
                m3=-(m1+m2)
                if(abs(m3) > l3)cycle
                wigner(m1,m2,m3)=threej(l1,l2,l3,m1,m2,m3)
            enddo
        enddo

    end function wigner

    !****************************************************************************
    !** Function calculating the Wigner 3-j coefficient                        **
    !****************************************************************************
    function threej(j1,j2,j3,m1,m2,m3)
        integer, intent(in) :: j1,j2,j3,m1,m2,m3
        real(PR) :: threej

        integer :: t1,t2,t3,t4,t5
        integer :: t,tmin,tmax

        if(j3 > j1 + j2 .or. j3 < abs(j1 - j2))then
            write(*,*)'j3 is out of bounds.'
            return
        endif

        if(abs(m1) > j1)then
            write(*,*)'m1 is out of bounds.'
            return
        endif

        if(abs(m2) > j2)then
            write(*,*)'m2 is out of bounds.'
            return
        endif

        if(abs(m3) > j3)then
            write(*,*)'m3 is out of bounds.'
            return
        endif

        t1 = j2 - m1 - j3
        t2 = j1 + m2 - j3
        t3 = j1 + j2 - j3
        t4 = j1 - m1
        t5 = j2 + m2

        tmin = max( 0, max( t1, t2 ) )
        tmax = min( t3, min( t4, t5 ) )

        threej = 0.0

        do t = tmin,tmax
            threej = threej + (-1)**t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) &
                * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )
        enddo

        threej = threej * (-1)**(j1-j2-m3) &
            * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1) &
            * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) )
    end function threej

    !*******************************************************************************
    !** Function calculating the factorial                                        **
    !*******************************************************************************
    function factorial(n)
        integer, intent(in) :: n
        real(PR) :: factorial

        integer, save :: ntop=0
        integer, parameter :: nmax=33
        real(PR), dimension(nmax), save :: a
        integer :: j

        a(1)=1.0
        if(n <= ntop)then
            factorial=a(n+1)
        elseif(n < nmax)then
            do j=ntop+1,n
                a(j+1)=j*a(j)
            enddo
            ntop=n
            factorial=a(n+1)
        else
            stop 'argument in factorial large'
        endif
    end function factorial
    !*******************************************************************************
end module spherical_harmonics
