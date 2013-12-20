SUBROUTINE ellipse(fipsi,nres,elldih,resnam,tekrok,asseq,stralf,minres,ncach,nch,namcor,ellfipsi)
USE tablice
IMPLICIT none
	
real(8)::ellfipsi(15000,2),fipsi(15000,3),t1(7),t2(7),dist,cost,fi,psi,fiel,&
psiel,fielmin,psielmin,q,t,tp,y,xmindist
integer ncach(0:50),nres,elldih,stralf,minres,nch,kp,m,tekrok,tdeg,tpdeg
character*1 resnam(15000)*3,asseq(15000),namcor*50
		
kp=index(namcor,' ')
	
!====================================================================
!..ellipse creation
!.......parametr q elipsy
	q=45.0*(PI/180)

	do 22 m=1,nch
	  if((ncach(m)-ncach(m-1)).lt.minres)goto 22

	if(elldih.eq.1)then
        open(17,file=namcor(1:(kp-1))//"."//sfx(m)//'.elldih',status='unknown')
	elseif(elldih.ne.0)then
	write(16,850)'elldih'
	stop
	endif

	if(fipsi(ncach(m-1)+1,2).le.0)then
	tdeg=180
	else
	tdeg=0
	endif
	t=tdeg*(PI/180)
	psiel=-127.28*SIN(q)*COS(t)-84*COS(q)*SIN(t)
	ellfipsi(1,2)=psiel
	if(elldih.eq.1)then
	write(17,810)ncach(m-1)+1,resnam(ncach(m-1)+1),psiel
	endif

	do 20 i=ncach(m-1)+2,ncach(m)-1	

	fi=fipsi(i,1)
	psi=fipsi(i,2)

	if(fi.le.0)then
		if(psi.le.0)then
		tpdeg=90
		elseif(psi.gt.0)then
		tpdeg=180
		else
		write(*,*)' '
		write(*,*)'warning!  wrong psi value'
		write(*,*)' '
		endif
	elseif(fi.gt.0)then
		if(psi.le.0)then
                tpdeg=0
                elseif(psi.gt.0)then
                tpdeg=270
                else
		write(*,*)' '
                write(*,*)'warning!  wrong psi value'
		write(*,*)' '
                endif
	else
	write(*,*)' '
	write(*,*)'warning!  wrong fi value'
	write(*,*)' '
	endif
	
	tpdeg=tpdeg-60
	tp=tpdeg*(PI/180)
	fiel=127.28*COS(q)*COS(tp)-84*SIN(q)*SIN(tp)
	psiel=-127.28*SIN(q)*COS(tp)-84*COS(q)*SIN(tp)
	xmindist=SQRT((fiel-fi)**2+(psiel-psi)**2)
	
DO tdeg=tpdeg,tpdeg+120,tekrok
	t=tdeg*(PI/180)
	fiel=127.28*COS(q)*COS(t)-84*SIN(q)*SIN(t)
	psiel=-127.28*SIN(q)*COS(t)-84*COS(q)*SIN(t)
	dist=SQRT(((fiel-fi)**2)+((psiel-psi)**2))
	if(dist.le.xmindist)then
	xmindist=dist
	fielmin=fiel
	psielmin=psiel
	endif 
END DO	 	

16	ellfipsi(i,1)=fielmin
	ellfipsi(i,2)=psielmin

	if(elldih.eq.1)then
	write(17,811)i,resnam(i),fielmin,psielmin
	endif

20	continue

	if(fipsi(ncach(m),1).le.0)then
        tdeg=180
        else
        tdeg=0
        endif
	t=tdeg*(PI/180)
	fiel=127.28*COS(q)*COS(t)-84*SIN(q)*SIN(t)
	ellfipsi(nres,1)=fiel

	if(elldih.eq.1)then
		write(17,812)nres,resnam(nres),fiel
      elseif(elldih.ne.0)then
          write(16,850)'elldih'
          stop
      endif

	close(17)
22	continue

!======================================================================
!......structural alfabet assignment - basing on previously established
!......alfabet definition
if(stralf.eq.1)then
  do i=1,7
	t1(i)=lettersdef(i,1)*(PI/180)
    t2(i)=lettersdef(i,2)*(PI/180)
  end do

  do m=1,nch
	if((ncach(m)-ncach(m-1)).lt.minres)cycle
	do k=ncach(m-1)+2,ncach(m)-1	
	  fi=ellfipsi(k,1)
	  psi=ellfipsi(k,2)
	  !	write(*,*)psi,fi
	  cost=-((fi-psi)/(2*127.28*cos(q)))
	  t=PI-ACOS(cost)
	  y=(psi-fi)/2
	  if(psi.gt.y)t=2*PI-t
	  do i=1,7
		if(t.ge.t1(i).AND.t.lt.t2(i))then
		asseq(k)=AS(i)
		!	 write(*,*)AS(i)
		endif
	  end do	
	end do
  end do

elseif(stralf.ne.0)then
  write(16,850)'stralf'
  stop
endif

asseq(1)=asseq(2)
asseq(nres)=asseq(nres-1)

810	format(i7,1x,a3,1x,8x,f8.3)
811	format(i7,1x,a3,1x,2f8.3)
812	format(i7,1x,a3,1x,f8.3)
813	format(i7,1x,i7,1x,a3,f8.3)
850	format('ERROR: Wrong value of parameter: ',a8,/&
       'restart with proper value set')

	
END
