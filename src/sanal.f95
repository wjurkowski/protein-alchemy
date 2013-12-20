SUBROUTINE sanal(atomy,fatom,bocz,fbocz,nres,fnres,fnat,ncach,nch,nseq,namcor,fnamcor,dccaprof,cmap,cmaptype,cmapcut,overlay,&
numf,bbaln,AAoverlay,nat,cmapmol,num_mol,resnam,resid,chainid,mostki,fit,capri,recs,rece,ligs,lige)
USE tablice
IMPLICIT none

REAL(8)::atomy(15000,30,4),fatom(15000,30,4),DcCadiff,dist,sumax,sumay,sumaz,x1,y1,z1,x2,y2,z2,rmsd,nofitrmsd,rmsdback,cdist,d,&
nofitdist,ct,ssdat(100,4)
INTEGER::bocz(15000),fbocz(15000),ncach(0:50),dccaprof,nch,nseq,cmap,cmaptype,cmapcut,kpos,kk,l,m,mm,n,nn,j,nat,blad,fnres,lss,&
fnat,nres,fm,bbaln,nr,fnr,overlay,kp1,kp2,numf,lt,ilt,nre,pier,osta,kres,AAoverlay,cmapmol,ctt,t1,t2,num_mol,wart(8),&
resid(15000),lsiar,mostki,fit,recs,rece,ligs,lige,capri
character namcor*50,fnamcor*50,data*8,czas*10,strefa*5,tab*1,resnam(15000)*3,chainid(15000)*3,ssdat2(100,2)*3
REAL(8), ALLOCATABLE :: X(:,:),Y(:,:),AArmsd(:),AAdist(:),nofitAArmsd(:),nofitAAdist(:),aXm(:,:),aYm(:,:),geomc(:,:),&
DcCa(:,:,:)
INTEGER, ALLOCATABLE :: fressrt(:),ressrt(:),fresend(:),resend(:),lkont(:),kontakt(:,:)

tab=char(9)
kpos=index(namcor,' ')
kp1=index(fnamcor,'.')
kp2=index(namcor,'.')

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2	
if(dccaprof.eq.1)then
  call DATE_AND_TIME(data,czas,strefa,wart)
  write(16,*)'calculates DcCalfa profiles...'
  !......Dc-Calfa profiles
  open(18,file=(namcor(1:(kpos-1))//'.DcCa'))
  ALLOCATE(geomc(nch,3),DcCa(nseq,nch,nres),stat=blad)
  if(blad.ne.0) write(*,*)'Problemy z tablicami:',blad
  write(18,'(A20,A9,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with prot_alchem',&
  ' on:',wart(1),'-',wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
  write(18,'(A35,A15,A15)')'Files analysed [PDB and chain ID]: ',namcor
  do m=1,nch
	l=0
	do k=ncach(m-1)+1,ncach(m)
	  do n=1,bocz(k)
		x1=atomy(k,n,2)
        y1=atomy(k,n,3)
        z1=atomy(k,n,4)
		sumax=sumax+x1
		sumay=sumay+y1
		sumaz=sumaz+z1
		l=l+1
	  end do
	end do
	geomc(m,1)=sumax/l
	geomc(m,2)=sumay/l
	geomc(m,3)=sumaz/l
  end do
  do m=1,nch
	do k=ncach(m-1)+1,ncach(m)
	  x1=atomy(k,2,2)
	  y1=atomy(k,2,3)
	  z1=atomy(k,2,4)
	  DcCa(nseq,m,k)=SQRT((x1-geomc(m,1))**2+(y1-geomc(m,2))**2+(z1-geomc(m,3))**2)
	end do
  end do

  if(nseq.eq.2)then
	do m=1,nch
	  do k=ncach(m-1)+1,ncach(m)
		DcCadiff=DcCa(1,m,k)-DcCa(2,m,k)
		write(18,800)m,k,DcCa(nseq,m,k),DcCadiff
	  end do
	end do
  elseif(nseq.eq.1)then
	do m=1,nch
	  do k=ncach(m-1)+1,ncach(m)
		write(18,801)m,k,DcCa(nseq,m,k)
	  end do
	end do
  endif
  close(18)
  write(16,*)'...profiles calculated'	
elseif(dccaprof.ne.0)then
  write(*,850)'dccaprof'
  stop
endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
if(cmap.eq.1)then
  write(16,*)'calculates contact maps...'	
  call DATE_AND_TIME(data,czas,strefa,wart)
  !......contact map
  ALLOCATE(lkont(nres),kontakt(nres,100),stat=blad)
  if(blad.ne.0) write(*,*)'Problemy z tablicami:',blad
  if(cmapmol.eq.1)then
	do m=1,nch
	  open(18,file=(namcor(1:(kpos-1))//sfx(m)//'.cmap.o'))
	  write(18,'(A20,A9,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with prot_alchem',&
	  ' on:',wart(1),'-',wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
	  write(18,'(A35,A15,A15)')'Files analysed [PDB and chain ID]: ',fnamcor,namcor
	  if(cmaptype.eq.1)then	
		do k=ncach(m-1)+1,ncach(m)-2
		  x1=atomy(k,2,2)
	      y1=atomy(k,2,3)
	      z1=atomy(k,2,4)
		  do kk=k+1,ncach(m)-1
			x2=atomy(kk,2,2)
	        y2=atomy(kk,2,3)
	        z2=atomy(kk,2,4)
			dist=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
			if(dist.le.cmapcut)then
			  write(18,802)k,kk
			else
			  write(18,802)k,0
			endif
		  end do	
		end do
	  elseif(cmaptype.eq.2)then
		do k=ncach(m-1)+1,ncach(m)-2
		  lkont(k)=0
		  do n=1,bocz(k)
			x1=atomy(k,n,2)
			y1=atomy(k,n,3)
			z1=atomy(k,n,4)
			do kk=k+1,ncach(m)-1
			  do nn=1,bocz(kk)
				x2=atomy(kk,nn,2)
				y2=atomy(kk,nn,3)
				z2=atomy(kk,nn,4)
				dist=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
				if(dist.le.cmapcut)then
				  lkont(k)=lkont(k)+1
				  kontakt(k,lkont(k))=kk
				endif
			  end do 
			end do	
		  end do
		  write(18,802)k,(kontakt(k,j),j=1,lkont(k))
		end do
	  else
		write(*,850)'cmaptype'
		stop
	  endif 
	  close(18)
	end do

  elseif(cmapmol.eq.2)then
	if(numf.gt.1)then
	  open(18,file=(namcor(1:(kpos-1))//sfx(m)//'.cmap.o'))
	  write(18,'(A20,A9,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with prot_alchem',&
	  ' on:',wart(1),'-',wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
	  write(18,'(A35,A15,A15)')'Files analysed [PDB and chain ID]: ',fnamcor,namcor
	  if(cmaptype.eq.1)then	
		do k=1,fnres
		  x1=fatom(k,2,2)
		  y1=fatom(k,2,3)
	      z1=fatom(k,2,4)
			do kk=1,nres
			  x2=atomy(kk,2,2)
			  y2=atomy(kk,2,3)
			  z2=atomy(kk,2,4)
			  dist=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
			  ct=dist/(rvdw(1)+rvdw(1))
			  if(ct.lt.0.7)ctt=1
			  if(ct.lt.0.9.AND.ct.ge.0.7)ctt=2
			  if(ct.lt.1.2.AND.ct.ge.0.9)ctt=3
			  if(dist.le.cmapcut)then
				write(18,802)k,kk
			  else
			 	write(18,802)k,0
			  endif
			end do	
	 	end do
	  elseif(cmaptype.eq.2)then
		do k=1,fnres
		  lkont(k)=0
		  do n=1,fbocz(k)
			t1=fatom(k,n,1)
			x1=fatom(k,n,2)
	        y1=fatom(k,n,3)
	        z1=fatom(k,n,4)
			do kk=1,nres
			  do nn=1,bocz(kk)
				t2=atomy(kk,nn,1)
		 		x2=atomy(kk,nn,2)
	        	y2=atomy(kk,nn,3)
	        	z2=atomy(kk,nn,4)
				dist=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
				ct=dist/(rvdw(t1)+rvdw(t2))
				if(ct.lt.0.7)ctt=1
				if(ct.lt.0.9.AND.ct.ge.0.7)ctt=2
				if(ct.lt.1.2.AND.ct.ge.0.9)ctt=3
				if(dist.le.cmapcut)then
				  lkont(k)=lkont(k)+1
				  kontakt(k,lkont(k))=kk
				endif
			  end do 
 	 		end do	
		  end do
		  write(18,802)k,(kontakt(k,j),j=1,lkont(k))
		end do
	  else
		write(*,850)'cmaptype'
		stop
	  endif 
	  close(18)
	endif
	
  else
	write(*,850)'cmapmol'
	stop
  endif 	

  write(16,*)'...contact maps calculated'
elseif(cmap.ne.0)then
  write(*,850)'cmap'
  stop
endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
if(overlay.eq.1)then
  if(numf.gt.1)then
	kres=MIN(fnres,nres)
	write(16,*)'overlays and compares molecules...'	
	ALLOCATE(X(3,nat),Y(3,fnat),fressrt(fnres),ressrt(nres),fresend(fnres),resend(nres),stat=blad)
	if(AAoverlay.eq.1)ALLOCATE(AArmsd(kres),AAdist(kres),nofitAArmsd(kres),nofitAAdist(kres),stat=blad)
	if(blad.ne.0) write(*,*)'Problemy z tablicami:',blad
	if(AAoverlay.eq.1)open(18,file=(fnamcor(1:(kp1-1))//"-"//namcor(1:(kp2-1))//'.rmsd'),status='unknown')
	write(16,*)'calculates rmsd...'
	call DATE_AND_TIME(data,czas,strefa,wart)
	if(AAoverlay.eq.1.AND.numf.eq.1)then
		write(18,'(A20,A9,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with Siatami',&
		' v8.8 on:',wart(1),'-',wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
		write(18,'(A35,A15,A15)')'Files analysed [PDB and chain ID]: ',fnamcor,namcor
		write(18,'(a)')'Selected fragments of aligned sequences were superimposed separately'
		write(18,'(a)')'Rmsd calculated for given motif (total RMSD value) and for each AA'
	endif
!......whole molecules
	if(bbaln.eq.0)then	!all atoms
		m=0
		fm=0
		do k=1,fnres
		  fressrt(k)=fm+1
		  do n=1,fbocz(k)
			fm=fm+1
			Y(1,fm)=fatom(k,n,2)
			Y(2,fm)=fatom(k,n,3)
			Y(3,fm)=fatom(k,n,4)
	 	  end do
		  fresend(k)=fm	
		end do
		do k=1,nres
		  ressrt(k)=m+1
		  do n=1,bocz(k)
			m=m+1
			X(1,m)=atomy(k,n,2)
			X(2,m)=atomy(k,n,3)
			X(3,m)=atomy(k,n,4)
		  end do
		  resend(k)=m
		end do
	endif
	if(bbaln.eq.1)then	!backbone only
		m=0
		fm=0
		do k=1,fnres
		  fressrt(k)=fm+1	
	   	  do n=1,4
			fm=fm+1
			Y(1,fm)=fatom(k,n,2)
			Y(2,fm)=fatom(k,n,3)
			Y(3,fm)=fatom(k,n,4)	
		  end do
		  fresend(k)=fm
 		end do
		do k=1,nres
		  ressrt(k)=m+1
		  do n=1,4
			m=m+1
			X(1,m)=atomy(k,n,2)
			X(2,m)=atomy(k,n,3)
			X(3,m)=atomy(k,n,4)
		  end do
		  resend(k)=m
		end do
	endif
	if(kres.eq.fnres)lt=fm
	if(kres.eq.nres)lt=m
	nofitdist=cdist(Y,X,lt)
	nofitrmsd=rmsdback(Y,X,lt)
	write(16,*)'finds best fit for molecules in question...'
	if(fit.eq.1)then
                call fitting(X,Y,lt,rmsd,dist)
	        write(*,804)fnamcor,tab,namcor,tab,rmsd,dist
        elseif(fit.eq.0)then
	        write(*,804)fnamcor,tab,namcor,tab,nofitrmsd,tab,nofitdist
        endif
	write(16,*)'...molecules compared'

	if(AAoverlay.eq.1)then          !single residue analysis
	do m=1,kres
		pier=fressrt(m)
		osta=fresend(m)
		if(nres.lt.fnres)then
		  pier=ressrt(m)
		  osta=resend(m)
		endif
		k=osta-pier+1
		ALLOCATE(aXm(3,k),aYm(3,k),stat=blad)
		if(blad.ne.0) write(*,*)'Problemy z tablicami AAoverlay:',blad
		k=0
		do i=pier,osta
		  k=k+1
	      aXm(1,k)=X(1,i)
		  aXm(2,k)=X(2,i)
		  aXm(3,k)=X(3,i)
		  aYm(1,k)=Y(1,i)
	      aYm(2,k)=Y(2,i)
	      aYm(3,k)=Y(3,i)
		end do
	    nofitAArmsd(m)=rmsdback(aYm,aXm,k)
	    nofitAAdist(m)=cdist(aYm,aXm,k)
	    if(fit.eq.1)call fitting(aXm,aYm,k,rmsd,dist)
	    AArmsd(m)=rmsd
	    AAdist(m)=dist
	    write(18,815)m,tab,AArmsd(m),tab,nofitAArmsd(m),tab,AAdist(m),tab,nofitAAdist(m)
	end do
	close(18)
	endif
  
  elseif(numf.eq.1)then
	write(16,*)'first structure, waiting for the second one' 
  else
	write(*,850)'numf'
	stop
  endif
elseif(overlay.ne.0)then
	write(*,850)'overlay'
	stop
endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
if(mostki.eq.1)then !.....counting of SS-bridges
  write(16,*)'calculates number of SS-bridges...'
  l=0
  lss=0
  do k=1,nres
	t1=atomy(kk,nn,1)
	if(atomy(kk,nn,1).eq.4)then
	  l=l+1
	  ssdat(l,1)=atomy(k,2,2)
	  ssdat(l,2)=atomy(k,2,3)
	  ssdat(l,3)=atomy(k,2,4)
	  ssdat(l,4)=resid(k)
	  ssdat2(l,1)=resnam(k)
	  ssdat2(l,2)=chainid(k)
	end if
  end do
  lsiar=l
  write(16,817)"SS-bond",tab,"RES ID",tab,"RESIDUE",tab,"CHAIN ID",tab,"RES ID",tab,"RESIDUE",tab,"CHAIN ID",tab
  do l=1,lsiar-1
	do m=l+1,lsiar
		  x1=ssdat(l,1)
		  y1=ssdat(l,2)
		  z1=ssdat(l,3)
		  x2=ssdat(m,1)
		  y2=ssdat(m,2)
		  z2=ssdat(m,3)
		  d=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		  if(d.ge.SSr1.AND.d.le.SSr2)then
			lss=lss+1
			write(16,816)lss,tab,ssdat(l,4),tab,ssdat2(l,1),tab,ssdat2(l,2),tab,ssdat(m,4),tab,ssdat2(m,1),tab,ssdat2(m,2)
		  endif
	end do
  end do
  write(16,*)'...SS-bridges calculated'
  write(16,'(A15,A1,I3)')namcor,lss
elseif(mostki.ne.0)then
  write(*,850)'mostki'
  stop
endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!.....analyses structures
if(capri.eq.1)then
     if(numf.gt.1)then
	call capri_rms(atomy,fatom,bocz,fbocz,nres,fnres,fnat,namcor,fnamcor,bbaln,nat,resnam,resid,recs,rece,ligs,lige)
     elseif(numf.eq.1)then
	write(16,*)'first structure, waiting for the second one' 
     else
	write(*,850)'numf'
	stop
     endif
elseif(capri.ne.0)then
     write(*,850)'capri'
     stop
endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
800	format(I2,1x,I5,1x,f7.3,1x,f7.3)
801	format(I2,1x,I5,1x,f7.3)
802	format(I5,1x,100I5)
803	format(I5)
804	format(A15,A1,A15,A1,f7.3,A1,f7.3)	
815	format(I7,A1,4(F7.3,A1))
816	format(2(I7,A1),2(A3,A1),I7,A1,2(A3,A1))
817	format(7(A8,A1))
850   format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')

END SUBROUTINE
