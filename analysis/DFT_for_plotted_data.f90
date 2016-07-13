!this is going to obtain the fourier coefficients for ""plotted data"", fft_data.dat
!we obtain fourier expansion coefficients up to the 10th
program DFT_for_plotted
use functions_module
implicit none
double precision,dimension(0:9) :: c,s
double precision,dimension(1:2,0:10000):: PLOT
double precision,dimension(0:10000):: nu, P
integer :: i,ierrd,j,MIN,MAX
double precision :: d_nu,d_c,d_s
double precision:: nu_i,nu_f,per,nu_p,nu_m
!double precision,parameter :: L=2.25d0*25.4d0,c_l=299.792458d0 !c_l speed of light
double precision,parameter :: pi= 3.14159265359d0
open(10, file='fft_data.dat', status= 'old')
open(11, file='fourier_transformed.dat', status='unknown')
write(11,*) 'k,','c(k),','s(k)'

!first we read the plotted data and find the initial nu nu_i and the final nu nu_f
!the data to read is assumed as an array of nu|P
do i=0,10
read(10,*,iostat=ierrd) PLOT
nu=PLOT(1,0:10000)
P=PLOT(2,0:10000)
if (ierrd/=0) exit
end do

MIN=1 !the row number that has the first row of numbers, counting from zero!

do i=MIN,1000000
    if(nu(i)>0) then
    cycle
    else
    MAX=i
    nu_f=nu(i-1)
    exit
    end if
end do

! defining and splitting the period, convert to \int_{-L}^{L}
per=nu_f-nu_i
nu_p=per/2.d0
nu_m=(-1.d0)*nu_p
d_nu=(nu_f-nu_i)/dble(MAX-MIN) !assuming the data plots are in equal nu intervals

! start Fourier expansion
!c(j)=(2/per)\int_{nu_m}^{nu_p}f(nu)cos(2*j*pi*nu/per)d nu, Simpson's formula
!s(j)=(2/per)\int_{nu_m}^{nu_p}f(nu)sin(2*j*pi*nu/per)d nu, Simpson's formula
do j=0,9
c(j)=0.d0
s(j)=0.d0
nu=nu_m
do i=MIN,MAX-2
d_c=P(i)*cos(2.d0*j*pi*(nu(i)-per/2)/per)+4.d0*P(i+1)*cos(2.d0*j*pi*(nu(i+1)-per/2)/per)+P(i+2)*cos(2.d0*j*pi*(nu(i+2)-per/2)/per)
c(j)=c(j)+ d_c/6.d0 * d_nu/per
d_s=P(i)*sin(2.d0*j*pi*(nu(i)-per/2)/per)+4.d0*P(i+1)*sin(2.d0*j*pi*(nu(i+1)-per/2)/per)+P(i+2)*sin(2.d0*j*pi*(nu(i+2)-per/2)/per)
s(j)=s(j)+ d_s/6.d0 * d_nu/per
nu = nu + d_nu
end do
write(11,*) j,c(j),s(j)
end do
end program DFT_for_plotted

