SUBROUTINE fitting(X,Y,ltab,rmsd,dist)
IMPLICIT none

REAL(8)::X(3,ltab),X1(3,ltab),Y(3,ltab),Y1(3,ltab),XT(ltab,3),A(3,3),U(3,3),VT(3,3),S(3),WORK(60),&
Z(3,3),Xp(3,ltab),xcent2,ycent2,zcent2,zcent1,xcent1,ycent1,sumx,sumy,sumz,suma,Xpt(3,ltab),rmsd,&
rmsdback,aXm(3,4),aYm(3,4),calfX(3,1),calfY(3,1),cdist,dist
INTEGER::k,l,m,n,LDVT,INFO,LDU,LDA,LWORK,ltab,nr,i,nre
character*1 JOBU,JOBVT

!...geometric center of first motive		
	sumx=0
	sumy=0
	sumz=0
	do 9 m=1,ltab	
        sumx=sumx+Y(1,m)
        sumy=sumy+Y(2,m)
        sumz=sumz+Y(3,m)
9	continue
	xcent1=sumx/ltab     
	ycent1=sumy/ltab
	zcent1=sumz/ltab

	do 10 m=1,ltab
	Y1(1,m)=Y(1,m)-xcent1
    	Y1(2,m)=Y(2,m)-ycent1
    	Y1(3,m)=Y(3,m)-zcent1		
10	continue		
	
	
!...geometric center of consecutive motive	
	sumx=0
	sumy=0
	sumz=0
		do 15 m=1,ltab	
        sumx=sumx+X(1,m)
        sumy=sumy+X(2,m)
        sumz=sumz+X(3,m)
15		continue
	xcent2=sumx/ltab
	ycent2=sumy/ltab
	zcent2=sumz/ltab
	

!...moves to the origin
		
	do 16 m=1,ltab
        X1(1,m)=X(1,m)-xcent2
        X1(2,m)=X(2,m)-ycent2
        X1(3,m)=X(3,m)-zcent2
16	continue					

!...finds A from  Y=A*X by transposing X 
!...Y*XT=A

	do 18 m=1,ltab
	XT(m,1)=X1(1,m)
	XT(m,2)=X1(2,m)
	XT(m,3)=X1(3,m)
18	continue

!...multiplies Y and XT	
	do 22 k=1,3
	do 21 m=1,3
	suma=0
	 do 20 n=1,ltab
	suma=suma+Y(m,n)*XT(n,k)
20	 continue
	A(m,k)=suma
21	continue
22	continue
!...calculates SVD of A
	JOBU='A'
	JOBVT='A'
	M=3
	N=3
	LDA=3
	LDU=3
	LDVT=3
	LWORK=60
	call DGESVD(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)
!	write(16,*)'INFO1: ',INFO

!..calculates U*VT
	do 40 l=1,3
         do 39 m=1,3
	 suma=0 
	  do 38 n=1,3
        suma=suma+U(m,n)*VT(n,l)
38	  continue
	Z(m,l)=suma	
39       continue
40      continue
	
!...multipies Z and X giving X1: the best molecular fit
	do 45 m=1,ltab
	 do 44 l=1,3
	 suma=0
	  do 43	n=1,3
	suma=suma+Z(l,n)*X1(n,m)
43	  continue	
	Xp(l,m)=suma	 
44	 continue
45	continue
	
	do 48 m=1,ltab	 
	Xpt(1,m)=Xp(1,m)+xcent2
	Xpt(2,m)=Xp(2,m)+ycent2
	Xpt(3,m)=Xp(3,m)+zcent2	 
48	continue

!#################END OF FIT
!...calculates rmsd
	dist=cdist(Y1,Xp,ltab)
	rmsd=rmsdback(Y1,Xp,ltab)
	!dist=cdist(Y1,Xp,ltab)	
!	write(*,*)dist,"WWWW",ltab

end subroutine

