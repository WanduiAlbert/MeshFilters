!first you define the function that you want to FT
module functions_module
implicit none
contains
double precision function f(nu)
double precision,parameter :: T=0.9763d0, R=0.0234d0, A=0.0003d0, L=2.25d0*25.4d0, c=299.792458d0
double precision,parameter :: pi= 3.14159265359d0
double precision:: nu_fsr, Fi, delta
double precision,intent(in) :: nu
nu_fsr = c/(2.d0*L)
Fi= 4.d0*R/(1.d0-R)**2
delta = 2.d0 * pi * nu/nu_fsr
f = (T/(1.d0-R))**2*1.d0/(1.d0+Fi*sin(delta/2.d0)**2)
end function
end module

!beginning of the program. we obtain fourier expansion coefficients up to the 10th
program DFT
use functions_module
implicit none
double precision,dimension(0:9)::c,s
integer :: i,j,N
double precision :: nu,d_nu,d_c,d_s
double precision:: nu_i,nu_f,per,nu_p,nu_m,nu_fsr
double precision,parameter :: L=2.25d0*25.4d0,c_l=299.792458d0
double precision, parameter :: pi= 3.14159265359d0
open(11, file='fourier_ideal5.dat', status='unknown')
write(11,*) 'k,','c(k),','s(k)'

! defining and splitting the period, convert to \int_{-L}^{L}
nu_fsr = c_l/(2.d0*L)
nu_i=141.63232647d0
nu_f=nu_i+nu_fsr
per=nu_f-nu_i
nu_p=per/2.d0
nu_m=(-1.d0)*nu_p
N= 10000 !number of splits
d_nu=(nu_f-nu_i)/dble(N)

! start Fourier expansion
!c(j)=(2/per)\int_{nu_m}^{nu_p}f(nu)cos(2*j*pi*nu/per)d nu, Simpson's formula
!s(j)=(2/per)\int_{nu_m}^{nu_p}f(nu)sin(2*j*pi*nu/per)d nu, Simpson's formula
do j=0,9
c(j)=0.d0
s(j)=0.d0
nu=nu_m
    do i=1,N
    d_c=f(nu)*cos(2.d0*j*pi*nu/per)+4.d0*f(nu+d_nu)*cos(2.d0*j*pi*(nu+d_nu)/per)+f(nu+2.d0*d_nu)*cos(2.d0*j*pi*(nu+2.d0*d_nu)/per)
    c(j)=c(j)+ d_c/6.d0 * d_nu/per
    d_s=f(nu)*sin(2.d0*j*pi*nu/per)+4.d0*f(nu+d_nu)*sin(2.d0*j*pi*(nu+d_nu)/per)+f(nu+2.d0*d_nu) *sin(2.d0*j*pi*(nu+2.d0*d_nu)/per)
    s(j)=s(j)+ d_s/6.d0 * d_nu/per
    nu = nu + d_nu
    end do
write(11,*) j,c(j),s(j)
end do
end program DFT

