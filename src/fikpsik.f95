!===============================================================================
SUBROUTINE fikpsik(fipsi,minres,ncach,nch,atomy,resid,resnam,namcor,dihsave,dih_outliers)
USE tablice
IMPLICIT none
!***********  program do obliczania katow dwusciennych z pdb ********************
!........ Wiktor Jurkowski

character*15 resnam(15000)*3,namcor*50,kwas*1,tab*1
logical nofi,nopsi,savetors,noomega
real(8)::fipsi(15000,3),atomy(15000,30,4),coordCA(3,0:2),coordN(3,0:2),coordC(3,0:2),fi1(4,4),&
fi2(4,4),psi1(4,4),psi2(4,4),ome1(4,4),ome2(4,4),fi,psi,omega  
integer ncach(0:50),minres,nch,dihsave,kpos,m,lp,dih_outliers,resid(15000)

tab=char(9)

savetors=.FALSE.
	fi=0.0
	psi=0.0
	omega=0.0
	kpos=index(namcor,'.')
	if(dihsave.eq.1)then
	  if(dih_outliers.eq.1)then
	  open(18,file=("nonstandard_phi.dat"),access='APPEND',status='unknown',action='readwrite',err=100)
	  open(19,file=("nonstandard_psi.dat"),access='APPEND',status='unknown',action='readwrite',err=100)
	  open(20,file=("nonstandard_omega.dat"),access='APPEND',status='unknown',action='readwrite',err=100)
	  endif
	endif
	do 60 m=1,nch
	  if((ncach(m)-ncach(m-1)).lt.minres)goto 60
	if(dihsave.eq.1)then
	  write(*,*)namcor(1:(kpos-1))//"."//sfx(m)//".dih"
	  open(17,file=(namcor(1:(kpos-1))//"."//sfx(m)//".dih"),status='unknown',err=100)
	  savetors=.TRUE.
	elseif(dihsave.ne.0)then
	write(*,850)'dihsave'
	stop
	endif

	 do 58 k=ncach(m-1)+1,ncach(m)
	coordN(1,1)=atomy(k,1,2)
	coordN(2,1)=atomy(k,1,3)
	coordN(3,1)=atomy(k,1,4)

	coordCA(1,1)=atomy(k,2,2)
	coordCA(2,1)=atomy(k,2,3)
	coordCA(3,1)=atomy(k,2,4)

	coordC(1,1)=atomy(k,3,2)
	coordC(2,1)=atomy(k,3,3)
	coordC(3,1)=atomy(k,3,4)
		
	coordN(1,2)=atomy(k+1,1,2)
	coordN(2,2)=atomy(k+1,1,3)
	coordN(3,2)=atomy(k+1,1,4)
	
	coordCA(1,2)=atomy(k+1,2,2)
	coordCA(2,2)=atomy(k+1,2,3)
	coordCA(3,2)=atomy(k+1,2,4)

	
!..  obliczanie fi i psi
	if(k.eq.ncach(m-1)+1)then
	  nofi=.TRUE.
	  nopsi=.FALSE.
	  noomega=.FALSE.
	endif  
	if(k.eq.ncach(m))then 
		nopsi=.TRUE.
		noomega=.TRUE.
	endif
	
	do 32 i=1,3
	fi1(1,i)=1
	fi2(1,i)=1
	psi1(1,i)=1
	psi2(1,i)=1
	ome1(1,i)=1
	ome2(1,i)=1
32	continue	

	do 34 i=1,4
        fi1(i,4)=1
        fi2(i,4)=1
        psi1(i,4)=1
        psi2(i,4)=1
		ome1(i,4)=1
        ome2(i,4)=1
34	continue
 
	if(nofi)goto 37
	do 36 i=1,3
	fi1(2,i)=coordC(i,0)
	fi1(3,i)=coordN(i,1)	
	fi1(4,i)=coordCA(i,1)
	fi2(2,i)=coordN(i,1)
	fi2(3,i)=coordCA(i,1)
	fi2(4,i)=coordC(i,1)
36	continue
	call winkiel(PI,fi1,fi2,fi)
	
	if(nopsi)goto 39 
37	do 38 i=1,3
	psi1(2,i)=coordN(i,1)
	psi1(3,i)=coordCA(i,1)
	psi1(4,i)=coordC(i,1)
	psi2(2,i)=coordCA(i,1)
	psi2(3,i)=coordC(i,1)
	psi2(4,i)=coordN(i,2)
38	continue
	call winkiel(PI,psi1,psi2,psi)

	do i=1,3
		ome1(2,i)=coordCA(i,1)
		ome1(3,i)=coordC(i,1)
		ome1(4,i)=coordN(i,2)
		ome2(2,i)=coordC(i,1)
		ome2(3,i)=coordN(i,2)
		ome2(4,i)=coordCA(i,2)
	end do 
	call winkiel(PI,ome1,ome2,omega)

!.. writing table with fi, psi values, for dihsave=1 saves also the values to file	
39	 do 15 i=1,20
	 if(resnam(k).eq.aminok(i))then
	 kwas=AA(i)
	 fipsi(k,3)=i
	 	if(nofi)then
		  fipsi(k,2)=psi
		  if(savetors)write(17,808)resid(k),tab,kwas,tab,tab,psi,tab,omega
		  if(psi.gt.psicut(1).and.psi.lt.psicut(2))write(19,811)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,tab,tab,&
psi,tab,omega
		  if(omega.gt.omegacut(1).and.omega.lt.omegacut(2))write(20,811)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,tab,&
tab,psi,tab,omega
	  	elseif(nopsi)then
	  	  fipsi(k,1)=fi
	  	  if(savetors)write(17,807)resid(k),tab,kwas,tab,fi
		  if(fi.gt.phicut(1).and.fi.lt.phicut(2))write(18,812)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,fi
	  	  lp=0
	  	else
	  	  fipsi(k,1)=fi
	  	  fipsi(k,2)=psi
	  	  if(savetors)write(17,809)resid(k),tab,kwas,tab,fi,tab,psi,tab,omega
		  if(fi.gt.phicut(1).and.fi.lt.phicut(2))write(18,810)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,fi,tab,psi,tab,&
		  omega
		  if(psi.gt.psicut(1).and.psi.lt.psicut(2))write(19,810)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,fi,tab,psi,&
tab,omega
		  if(omega.gt.omegacut(1).and.omega.lt.omegacut(2))write(20,810)namcor(1:(kpos-1)),tab,resid(k),tab,kwas,tab,fi,&
tab,psi,tab,omega
	  	endif
	 endif
15	 continue	

807	format(i7,1a,a3,1a,f8.3)
808	format(i7,1a,a3,1a,10x,2(1a,f8.3))
809	format(i7,1a,a3,3(1a,f8.3))
810	format(a12,1a,i7,1a,a3,3(1a,f8.3))
811	format(a12,1a,i7,1a,a3,1a,1a,3x,2(1a,f8.3))
812	format(a12,1a,i7,1a,a3,1a,f8.3)
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
	nofi=.FALSE.
!	nopsi=.FALSE.
!	noomega=.FALSE.

	do 42 i=1,3	
	coordN(i,0)=coordN(i,1)
	coordCA(i,0)=coordCA(i,1)
	coordC(i,0)=coordC(i,1)
!	coordN(i,1)=coordN(i,2)
!	coordCA(i,1)=coordCA(i,2)
!	coordC(i,1)=coordC(i,2)
42	continue
		
!...next residue
58	 continue
	close(17)
60	continue
	close(18)
	close(19)
	close(20)

70	write(16,*)'done.... '
goto 101
100 write(*,*)'ERROR cant opene the file'
101      end

