SUBROUTINE blocker(fnamcor,namcor,atomy,fatom,ncach,nch,seq_motif2,fnat,numf,resid,fresid,nmot,fnres,nres,alnall,aaseq,faaseq,&
ignmtf,AAoverlay)
USE tablice
IMPLICIT none

REAL(8)::atomy(15000,30,4),fatom(15000,30,4),rmsd,ran,dist,nofitdist,nofitrmsd,cdist,rmsdback
INTEGER::ncach(0:50),j,l,m,fm,n,kp1,kp2,lt,nch,nr,fnat,num_mol,blad,w1,numf,nmot,AAoverlay,lltab,resid(15000),fresid(15000),&
fnres,nres,alnall,ktab,fnr,wart(8),dn,pier,osta,icutoff(20,10,2),aanumdif,ignmtf
character fnamcor*50,namcor*50,data*8,czas*10,strefa*5,aaseq(15000)*1,faaseq(15000)*1,aaseq2(5000,2)*1,tab*1
REAL(8), ALLOCATABLE :: X(:,:),Y(:,:),aXm(:,:),aYm(:,:),AArmsd(:),AAdist(:),nofitAArmsd(:),nofitAAdist(:)
INTEGER, ALLOCATABLE :: seq_motif2(:,:,:),resid2(:,:),fressrt(:),ressrt(:),fresend(:),resend(:)
LOGICAL gaps

!TYPE wynik
!	character fnamcor*50,namcor*50
!	real :: suff=7
!END TYPE
!TYPE(wynik) :: wynikowy

!suff=RAND(7)
!write(*,*)wynikowy%fnamcor,wynikowy%namcor,wynikowy%suff
	
kp1=index(namcor,' ')
write(16,*)'inter-molecular motifs alignment'
!paste sequence motifs
ALLOCATE(X(3,fnat),Y(3,fnat),seq_motif2(nmot,20,2),resid2(nres,2),fressrt(fnres),ressrt(nres),fresend(fnres),&
resend(nres),stat=blad)
if(blad.ne.0) write(*,*)'Problemy z tablicami X or Y:',blad
kp1=index(fnamcor,'.')
kp2=index(namcor,'.')
if(AAoverlay.eq.1)then
  ALLOCATE(AArmsd(fnres),AAdist(fnres),nofitAArmsd(fnres),nofitAAdist(fnres),stat=blad)
  if(blad.ne.0) write(*,*)'Problemy z tablicami AAoverlay:',blad 
endif
if(AAoverlay.eq.1)then
  open(18,file=(fnamcor(1:(kp1-1))//"-"//namcor(1:(kp2-1))//"."//sfx(numf)//'.rmsd'),status='unknown')
  call DATE_AND_TIME(data,czas,strefa,wart)
  if(alnall.eq.0)then
	write(18,'(A20,A9,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with Siatami',&
	' v8.8 on:',wart(1),'-',wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
	write(18,'(A35,A15,A15)')'Files analysed [PDB and chain ID]: ',fnamcor,namcor
	write(18,'(a)')'Selected fragments of aligned sequences were superimposed separately'
	write(18,'(a)')'Rmsd calculated for given motif (total RMSD value) and for each AA'
  elseif(alnall.eq.1)then
	write(18,'(A30,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Created with Siatami v8.8 on:',wart(1),'-',&
	wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
	write(18,*)'Selected fragments of aligned sequences were merged and superimposed'
	write(18,*)'Rmsd calculated for the resulted block (total RMSD value) and for each AA separately'
  endif
endif

write(16,*)'calculates rmsd...'
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!checks integrity of chain and detects gaps in structure
!first structure
gaps=.FALSE.
dn=0
do i=2,fnres
  aanumdif=fresid(i)-fresid(i-1)
  if(aanumdif.ne.1)then
	 gaps=.TRUE.
	 dn=dn+1
	 write(*,*)'WARNING:structure not present in sequence #: ',numf,'between AAs:',resid(i),resid(i-1)
	 write(16,'(A)')'WARNING:structure not present in sequence #: ',numf,'between AAs:',resid(i),resid(i-1)
	 icutoff(1,dn,1)=fresid(i-1)+1
	 icutoff(1,dn,2)=fresid(i)-1
  endif
end do
if(gaps)then
  do k=1,dn
	write(*,*)'GAP #: ',k, 'gap definition: ', icutoff(1,k,1),icutoff(1,k,2)
	do j=1,nmot
	  if(icutoff(1,k,1).le.seq_motif2(j,1,2).AND.icutoff(1,k,2).ge.seq_motif2(j,1,1))then
		write(*,*)'AFFECTED: motif #: ',j,'motif definition: ',	seq_motif2(j,1,1),seq_motif2(j,1,2)
	  endif	
	end do
  end do
endif
!next structure
gaps=.FALSE.
dn=0
do i=2,nres
  aanumdif=resid(i)-resid(i-1)
  if(aanumdif.ne.1)then
	gaps=.TRUE.
	dn=dn+1
	write(*,*)'WARNING:structure not present in sequence #: ',numf,'between AAs:',resid(i),resid(i-1)
	  write(16,'(A)')'WARNING:structure not present in sequence #: ',numf,'between AAs:',resid(i),resid(i-1)
	  icutoff(numf,dn,1)=resid(i-1)+1
	  icutoff(numf,dn,2)=resid(i)-1
  endif
end do
if(gaps)then
  do k=1,dn
	write(*,*)'GAP #: ',k, 'gap definition: ', icutoff(numf,dn,1),icutoff(numf,dn,2)
	do j=1,nmot
	  if(icutoff(numf,k,1).le.seq_motif2(j,numf,2).AND.icutoff(numf,k,2).ge.seq_motif2(j,numf,1))then
		write(*,*)'AFFECTED: part of motif #: ',j,'motif definition: ',	seq_motif2(j,numf,1),seq_motif2(j,numf,2)
	  endif
	end do
  end do
endif
  
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
m=0
fm=0
fnr=0
nr=0
lltab=0
do j=1+ignmtf,nmot
  do k=seq_motif2(j,1,1),seq_motif2(j,1,2)
	fressrt(k)=fm+1
	fnr=fnr+1
	do n=1,4
	  fm=fm+1
	  Y(1,fm)=fatom(k,n,2)
	  Y(2,fm)=fatom(k,n,3)
	  Y(3,fm)=fatom(k,n,4)
	  resid2(fnr,1)=k
	  aaseq2(fnr,1)=faaseq(k)
	end do    	
	fresend(k)=fm	
  end do
  do k=seq_motif2(j,numf,1),seq_motif2(j,numf,2)
	ressrt(k)=m+1
	nr=nr+1
	do n=1,4
	  m=m+1
	  X(1,m)=atomy(k,n,2)
	  X(2,m)=atomy(k,n,3)
	  X(3,m)=atomy(k,n,4)
	  resid2(nr,2)=k
	  aaseq2(nr,2)=aaseq(k)
	end do    	
	resend(k)=m	
  end do
			
  if(alnall.eq.0)then
	write(16,*)'finds best fit for molecules in question...'
	ktab=fm-lltab
	nofitdist=cdist(Y,X,ktab)
	nofitrmsd=rmsdback(Y,X,ktab)
	call fitting(X,Y,ktab,rmsd,dist)
	write(16,*)'...best fit found'
	write(16,*)'calculates rmsd...'
	write(*,816)j,tab,fnamcor,tab,namcor,tab,rmsd,tab,nofitrmsd,tab,dist,tab,nofitdist
	write(16,*)'...molecules compared'
	if(AAoverlay.eq.1)then          !single residue analysis
	  do m=1,fnr
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
	    call fitting(aXm,aYm,k,rmsd,dist)
	    AArmsd(m)=rmsd
	    AAdist(m)=dist
	    write(18,815)m,tab,AArmsd(m),tab,nofitAArmsd(m),tab,AAdist(m),tab,nofitAAdist(m)
	  end do
	endif
  endif
	lltab=fm		
end do 

if(alnall.eq.0)then
  close(18)
  write(16,*)'...rmsd calculated'
  write(16,*)'...alignment finished'
endif

!sprawdzenie dlugosci motywow
if(fm.ne.m)then
  write(*,*)'WARNING: different lenght of aligned motifs',fm,m
  write(16,*)'WARNING: different lenght of aligned motifs'
  lt=MIN(fm,m)
endif

!porownanie calkowitej dlugosci z iloscia AA w PDB
if(fnr.ne.fnres)then
  write(*,*)'WARNING: DELETION in sequence 1, different lenght of merged motifs and PDB',fnr,fnres
endif
if(nr.ne.nres)then
  write(*,*)'WARNING: DELETION in sequence ',numf,', different lenght of merged motifs and PDB',nr,nres
endif

!###########################################
if(alnall.eq.1)then
  write(16,*)'finds best fit for molecules in question...'
  nofitdist=cdist(Y,X,ktab)
  nofitrmsd=rmsdback(Y,X,ktab)
  call fitting(X,Y,lt,rmsd)
  write(16,*)'...best fit found'
  write(16,*)'calculates rmsd...'
  write(*,816)j,tab,fnamcor,tab,namcor,tab,rmsd,tab,nofitrmsd,tab,dist,tab,nofitdist
	write(16,*)'...molecules compared'
	if(AAoverlay.eq.1)then          !single residue analysis
	  do m=1,fnr
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
	    call fitting(aXm,aYm,k,rmsd,dist)
	    AArmsd(m)=rmsd
	    AAdist(m)=dist
	    write(18,815)m,tab,AArmsd(m),tab,nofitAArmsd(m),tab,AAdist(m),tab,nofitAAdist(m)
	  end do
	endif
  close(18)
  write(16,*)'...rmsd calculated'
  write(16,*)'...alignment finished'
endif

100	write(16,*)'...done'
814	format(8(a11))
815	format(I11,I11,A11,I11,A11,3(f11.6))
816	format(I3,A1,A15,A1,A15,A1,f7.3,A1,f7.3,A1,f7.3,A1,f7.3)
850	format('ERROR: Wrong value of parameter: ',a8,/&
       'restart with proper value set')
end
