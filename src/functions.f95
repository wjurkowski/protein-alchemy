!....various useful functions
!........rmsd all atoms

REAL(8) function rmsdall(A,B,bocz,ncach,nch)

REAL(8)::A(15000,30,4),B(15000,30,4),sq,x1,y1,z1,x2,y2,z2,rmsdall,sumsq
INTEGER::bocz(15000),ncach(0:50),k,l,m,n,nch
	l=0
	sumsq=0
	do 22 m=1,nch
		do 20 k=ncach(m-1)+1,ncach(m)
	 	 do 18 n=1,bocz(k)
	x1=A(k,n,2)
        y1=A(k,n,3)
        z1=A(k,n,4)
        x2=B(k,n,2)
        y2=B(k,n,3)
        z2=B(k,n,4)
	sq=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
	sumsq=sumsq+sq
	l=l+1	
18	 	 continue 
20 		continue
22	continue

	rmsdall=SQRT(sumsq/l)
	
end
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!........rmsd
REAL(8) function rmsdback(A,B,j)  
REAL(8)::A(3,j),B(3,j),x1,x2,y1,y2,z1,z2,sq,sumsq
sumsq=0
do k=1,j
  x1=A(1,k)
  y1=A(2,k)
  z1=A(3,k)
  x2=B(1,k)
  y2=B(2,k)
  z2=B(3,k)
  sq=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  sumsq=sumsq+sq
end do
rmsdback=SQRT(sumsq/j)
END
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.......cartesian
REAL(8) function dcart(x1,x2,y1,y2)
REAL(8) x1,x2,y1,y2
dcart=SQRT((x1-x2)**2+(y1-y2)**2)
END
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.......manhattan
REAL(8) function dmanh(x1,x2,y1,y2)
REAL(8) x1,x2,y1,y2
dmanh=ABS(x1-x2)+ABS(y1-y2)
END
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.......czebyszew
REAL(8) function dczeb(x1,x2,y1,y2)
REAL(8) x1,x2,y1,y2,a,b
a=ABS(x1-x2)
b=ABS(y1-y2)
dczeb=MAX(a,b)
END
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.......canberra
REAL(8) function dcanb(x1,x2,y1,y2)
REAL(8) x1,x2,y1,y2,a,b
a=ABS((x1-x2)/(x1+x2))
b=ABS((y1-y2)/(y1+y2))
dcanb=a+b
END
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

REAL(8) function cdist(A,B,j)
REAL(8)::A(3,j),B(3,j),x1,x2,y1,y2,z1,z2,cx1,cx2,cy1,cy2,cz1,cz2,sumx1,sumx2,sumy1,sumy2,sumz1,sumz2
!odleglosci od srodkow geometrycznych backbonu reszt i wiekszych fragmentów
sumx1=0
sumy1=0
sumz1=0
sumx2=0
sumy2=0
sumz2=0
do k=1,j
	x1=A(1,k)
	y1=A(2,k)
	z1=A(3,k)
	x2=B(1,k)
	y2=B(2,k)
	z2=B(3,k)
	sumx1=sumx1+x1
	sumy1=sumy1+y1
	sumz1=sumz1+z1
	sumx2=sumx2+x2
	sumy2=sumy2+y2
	sumz2=sumz2+z2
!	write(*,*)'eeeee',k,x1,y1,z1,x2,y2,z2,sumx1,sumx2,sumy1,sumy2,sumz1,sumz2
end do
	cx1=sumx1/j
	cx2=sumx2/j
	cy1=sumy1/j
	cy2=sumy2/j
	cz1=sumz1/j
	cz2=sumz2/j

	cdist=SQRT((cx1-cx2)**2+(cy1-cy2)**2+(cz1-cz2)**2)
!write(*,*)cdist,"AAAAAA"
END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!.......w coeficients

REAL(8) function wcoef(nf,lsek)
REAL(8) val,wcoef
INTEGER nf,lsek
	
	val=((nf+1)/((lsek/7.0)+1))
	wcoef=10*LOG10(val)
!	write(*,*)'ll',nf,lsek,val,wcoef
	end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!....geometric center for each chain separately
real function geomc(A,j)
REAL(8)::A(3,j),geomc(3),x1,y1,z1,sumx1,sumy1,sumz1
sumx1=0
sumy1=0
sumz1=0
do k=1,j-1
	x1=A(1,k)
	y1=A(2,k)
	z1=A(3,k)
	sumx1=sumx1+x1
	sumy1=sumy1+y1
	sumz1=sumz1+z1
end do
geomc(1)=sumx1/j
geomc(2)=sumy1/j
geomc(3)=sumz1/j
	
END
