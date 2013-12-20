SUBROUTINE capri_rms(atomy,fatom,bocz,fbocz,nres,fnres,fnat,namcor,fnamcor,bbaln,nat,resnam,resid,recs,rece,ligs,lige)
USE tablice
IMPLICIT none

REAL(8)::atomy(15000,30,4),fatom(15000,30,4),rmsd,nofitrmsd,rmsdback,suma,cent1(3),cent2(3)
INTEGER::bocz(15000),fbocz(15000),nseq,kpos,kk,l,m,mm,n,nn,j,nat,blad,fnres,fnat,nres,fm,bbaln,&
lt,kres,wart(8),resid(15000),recs,rece,ligs,lige
character namcor*50,fnamcor*50,data*8,czas*10,strefa*5,tab*1,resnam(15000)*3
REAL(8), ALLOCATABLE :: X(:,:),Y(:,:),X1(:,:),Y1(:,:),X2(:,:),Z(:,:),geomc(:,:)
tab=char(9)

kres=MIN(fnres,nres)
write(16,*)'overlays and compares molecules...'	
ALLOCATE(X(3,nat),Y(3,fnat),Z(3,3),stat=blad)
if(blad.ne.0) write(*,*)'Problems with allocation of X,Y,fressrt,ressrt,fresend,resend,Z',blad
write(16,*)'calculates rmsd...'
call DATE_AND_TIME(data,czas,strefa,wart)

!#############################
if(bbaln.eq.0)then	!all atoms
	m=0
	fm=0
	do k=recs,rece
	  do n=1,fbocz(k)
		fm=fm+1
		Y(1,fm)=fatom(k,n,2)
		Y(2,fm)=fatom(k,n,3)
		Y(3,fm)=fatom(k,n,4)
 	  end do
	end do
	do k=recs,rece
	  do n=1,bocz(k)
		m=m+1
		X(1,m)=atomy(k,n,2)
		X(2,m)=atomy(k,n,3)
		X(3,m)=atomy(k,n,4)
	  end do
	end do
endif
if(bbaln.eq.1)then	!backbone only
	m=0
	fm=0
	do k=recs,rece
   	  do n=1,4
		fm=fm+1
		Y(1,fm)=fatom(k,n,2)
		Y(2,fm)=fatom(k,n,3)
		Y(3,fm)=fatom(k,n,4)	
	  end do
	end do
	do k=recs,rece
		do n=1,4
			m=m+1
			X(1,m)=atomy(k,n,2)
			X(2,m)=atomy(k,n,3)
			X(3,m)=atomy(k,n,4)
		end do
	end do
endif
if(kres.eq.fnres)lt=fm
if(kres.eq.nres)lt=m

call fitting_capri(X,Y,lt,cent1,cent2,Z)

!##########################	
if(bbaln.eq.0)then	!all atoms
	m=0
	fm=0
	do k=ligs,lige
		do n=1,fbocz(k)
			fm=fm+1
!			Y(1,fm)=fatom(k,n,2)
!			Y(2,fm)=fatom(k,n,3)
!			Y(3,fm)=fatom(k,n,4)
	 	end do
        end do
	do k=ligs,lige
		do n=1,bocz(k)
			m=m+1
!			X(1,m)=atomy(k,n,2)
!			X(2,m)=atomy(k,n,3)
!			X(3,m)=atomy(k,n,4)
		end do
	end do
endif
if(bbaln.eq.1)then	!backbone only
	m=0
        fm=0
	do k=ligs,lige
	   	do n=1,4
			fm=fm+1
			Y(1,fm)=fatom(k,n,2)
			Y(2,fm)=fatom(k,n,3)
			Y(3,fm)=fatom(k,n,4)	
		end do
	end do
	do k=ligs,lige
		do n=1,4
			m=m+1
			X(1,m)=atomy(k,n,2)
			X(2,m)=atomy(k,n,3)
			X(3,m)=atomy(k,n,4)
		end do
	end do
endif
if(kres.eq.fnres)lt=fm
if(kres.eq.nres)lt=m
	
ALLOCATE(X1(3,lt),Y1(3,lt),X2(3,lt),stat=blad)
if(blad.ne.0) write(*,*)'Problems with allocation of X1 or Y1',blad

do m=1,lt
    X1(1,m)=X(1,m)-cent2(1)
    X1(2,m)=X(2,m)-cent2(2)
    X1(3,m)=X(3,m)-cent2(3)
end do				

do m=1,lt
	Y1(1,m)=Y(1,m)-cent1(1)
        Y1(2,m)=Y(2,m)-cent1(2)
    	Y1(3,m)=Y(3,m)-cent1(3)		
end do		

do m=1,lt
        do l=1,3
	        suma=0
	        do n=1,3
	                suma=suma+Z(l,n)*X1(n,m)
                end do	
	        X2(l,m)=suma	 
        end do
end do
	
rmsd=rmsdback(Y1,X2,lt)

write(*,804)fnamcor,tab,namcor,tab,rmsd
  
804	format(A15,A1,A15,A1,f7.3)

END SUBROUTINE
