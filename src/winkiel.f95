SUBROUTINE winkiel(PI,w1,w2,xkat)
IMPLICIT none
	
REAL(8)::w1(4,4),w2(4,4),wsp1(4),wsp2(4),a(3,3),b(3,3),c(3),d(3),PI,xkat,xkatrad,cokat,WM,WC,WD,summ,&
suml,summ1,summ2
INTEGER::i,j,k,ll,n
!..	rozwijanie wyznacznikow

	do 30 k=1,4
		ll=0
	n=1+k	
	   do 25 j=1,4
	if(j.eq.k)goto 25
		ll=ll+1 
	   do 20 i=1,3 
	a(i,ll)=w1((i+1),j)
        b(i,ll)=w2(i+1,j)
20	   continue
25	   continue
		
		do 26 i=1,3
26		continue		

WM=(a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3)-a(3,1)*a(2,2)*a(1,3)-&
a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3))
	wsp1(k)=WM*((-1)**n)
WM=(b(1,1)*b(2,2)*b(3,3)+b(2,1)*b(3,2)*b(1,3)+b(3,1)*b(1,2)*b(2,3)-b(3,1)*b(2,2)*b(1,3)-&
b(1,1)*b(3,2)*b(2,3)-b(2,1)*b(1,2)*b(3,3))
        wsp2(k)=WM*((-1)**n)
30	continue


	suml=0.0
	summ1=0.0
	summ2=0.0
	do 60 k=1,3
	suml=suml+(wsp1(k)*wsp2(k))
	summ1=summ1+(wsp1(k)**2)
	summ2=summ2+(wsp2(k)**2)
60	continue
	summ=SQRT(summ1)*SQRT(summ2)
	cokat=suml/summ

	if(cokat.lt.-1.0)cokat=-1.0
	if(cokat.gt.1.0)cokat=1.0

	xkatrad=ACOS(cokat)
	xkat=(xkatrad*180)/PI
!......................	ustalanie znaku kata...................................c

!	OP=(wsp1(1)*a(3,1)+wsp1(2)*a(3,2)+wsp1(3)*a(3,3)+wsp1(4))
!     &	/(SQRT(summ1))	

!..	iloczyn wektorowy	

	c(1)=wsp1(2)*wsp2(3)-wsp1(3)*wsp2(2)
	c(2)=wsp1(3)*wsp2(1)-wsp1(1)*wsp2(3)
	c(3)=wsp1(1)*wsp2(2)-wsp1(2)*wsp2(1)

!..	wektor przeciecia plaszczyzn

	d(1)=a(3,1)-a(2,1)
	d(2)=a(3,2)-a(2,2)
	d(3)=a(3,3)-a(2,3)

!..	sprawdzenie komplanarnosci
!..	gdy WC	ma znak inny niz WD to kat ustawiony jako ujemny 

WC=wsp1(1)*wsp2(2)*c(3)+wsp2(1)*c(2)*wsp1(3)+c(1)*wsp1(2)*wsp2(3)-wsp1(3)*wsp2(2)*c(1)-&
wsp2(3)*c(2)*wsp1(1)-c(3)*wsp1(2)*wsp2(1)

WD=wsp1(1)*wsp2(2)*d(3)+wsp2(1)*d(2)*wsp1(3)+d(1)*wsp1(2)*wsp2(3)-wsp1(3)*wsp2(2)*d(1)-&
wsp2(3)*d(2)*wsp1(1)-d(3)*wsp1(2)*wsp2(1)
	
	if((WC.lt.0.AND.WD.gt.0).OR.(WC.gt.0.AND.WD.lt.0))then
	xkat=-xkat
	endif
	


!	if(OP.gt.0)then
!	xkat=-xkat
!	endif

	end

