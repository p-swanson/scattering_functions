PROGRAM Sk
implicit none

complex :: im = (0,1)
real,external :: dist
integer :: i,j,time_step, n_atoms,box_len,n_row,n_col,k_val,m


real, dimension(225,3) :: array_in
complex,  dimension(50625,1) :: array_mid
!real,  dimension(225) :: array_mid
real, dimension(5000,1) :: k,array_out

box_len = 90

n_atoms = 225

time_step = 1000

DO k_val = 1,5000
k(k_val,1)= (k_val*0.001)
ENDDO

open(1,file='f_k_scatter_data.txt')
PRINT*,'open file'
DO n_row = 1,225
read(1,*) (array_in(n_row,n_col), n_col=1,3)
ENDDO
close(1)
PRINT *, 'Read array'

DO m = 1,5000
DO i = 1,n_atoms
DO j = 1,n_atoms
IF (j .NE. i) THEN
array_mid(j,1) = EXP(im*k(m,1)*dist(array_in(i,:),array_in(j,:),box_len))
ELSE
CONTINUE
ENDIF
ENDDO
ENDDO
array_out(m,1) = 1+SUM(REAL(array_mid))/n_atoms
ENDDO

PRINT *, 'Loop finished'
open(2,file='S_k_data_fort.txt')
DO i=1,time_step
write(2,*) k(i,1),array_out(i,1)
ENDDO
close(2)
END PROGRAM Sk

FUNCTION dist(input_array1,input_array2,box_len)
implicit none
real, dimension(1,3) ,intent(in) :: input_array1,input_array2
real :: dx,dy,dz,dist
integer, intent(in) :: box_len
real :: halfL



halfL = box_len/2

dx = input_array2(1,1)-input_array1(1,1)
dy = input_array2(1,2)-input_array1(1,2)
dz = input_array2(1,3)-input_array1(1,3)


dist = SQRT(dx**2+dy**2+dz**2)

END FUNCTION dist
