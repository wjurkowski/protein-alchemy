SUBROUTINE eceppanal(emesh,enertab)
USE tablice
IMPLICIT none

!..reads in ecepp energy maps with defined mesh

REAL(8)::enertab(20,360,360),emesh
INTEGER::numbc,numbr,n
character line*80,infile*30
numbr=360/emesh
numbc=360/emesh
do k=1,20
!	if(k.eq.15)numbc=2

  infile=('main_out.'//aminok(k))
  open(7,file=infile,status='old',err=98)
5     read(7,'(a)',end=100,err=99)line
  if(line(1:3).eq.'psi')then
	read(7,'(a)',end=100,err=99)line
	do n=1,numbc
        read(7,*,end=100,err=99)(enertab(k,i,n),i=1,numbr)
	end do
  else
    goto 5
  endif
  close(7)
end do

goto 100
98	write(*,*)'ERROR: cannot open file: ',infile
stop
99      write(*,*)'ERROR while reading input file: ',infile
stop
100     end

