SUBROUTINE readalign(inpalign,kwsize,resid,wind2,nseq,nres,seq_motif2,algn_typ,nmot,bcutoff,ecutoff,&
ignmtf)
USE tablice
IMPLICIT none

character line*90,lines(500)*82,inpalign*50,signstab(4,10000)*1
integer::nres(100),nre(100,10000),resid(15000),kp2,kwsize,nseq,j,l,m,n,ca,cb,cc,cd,ce,cf,cg,nf,&
npoz,licz,kp,nrows,blad,seq_motif(100,20,2),numrs(20),beg(20),kon(20),motl,algn_typ,nmot,nl,ns,nr,kl,a,b,&
bcutoff(20),ecutoff(20),seq_motif2(100,20,2),ignmtf
real(8)::wnp(10000),wind(10000),wnp2(100,10000),wind2(100,10000),w,suma,wcoef
logical nul,full

EXTERNAL wcoef
			
!===========================================================
!..reads in FASTA (Pearson) format
kp2=index(inpalign,' ')

if(algn_typ.eq.1)then
!structural alignments
  open(8,file=inpalign(1:(kp2-1))//'.wc',status='unknown',err=100)	
  licz=INT(kwsize/2)
  i=0
5	read(7,'(a)',end=42,err=99)line
  i=i+1
  kp=LEN(line)
  !	if(i.eq.1)lseq=kp-1
  do j=1,kp-1
	signstab(i,j)=line(j:j)
!	write(*,'(a)')'kkk',signstab(1,j)
  end do
  goto 5

42	nrows=i	
  !...checks number of sequences
  nseq=0
  do l=4,nrows
	if(signstab(l,32).ne.' ')nseq=nseq+1
	if(signstab(l,1).eq.' '.AND.nseq.ne.0)goto 50
  end do
50	npoz=0
  do n=1,nseq
	nres(n)=0
  end do
  
!...analyses the structural alignments
  do l=4,nrows,nseq+2	
	do k=32,81
	  npoz=npoz+1
	  ca=0
	  cb=0
	  cc=0
	  cd=0
	  ce=0
	  cf=0
	  cg=0
	  if(signstab(l,k).eq.' ')then
		npoz=npoz-1
		cycle
	  endif
	  n=0
	  do i=l,l+nseq-1
		n=n+1
		if(signstab(i,k).ne.'-')then
		  nres(n)=nres(n)+1
		  nre(n,npoz)=nres(n)
		elseif(signstab(i,k).eq.'-')then
		  nre(n,npoz)=0
		endif
		if(i.eq.l)write(*,*)nres(1)
		!	 write(*,*)'dddd',n,npoz,nres(n),nre(n,npoz)
		if(signstab(i,k).eq.'A')ca=ca+1
		if(signstab(i,k).eq.'B')cb=cb+1
		if(signstab(i,k).eq.'C')cc=cc+1
		if(signstab(i,k).eq.'D')cd=cd+1
		if(signstab(i,k).eq.'E')ce=ce+1
		if(signstab(i,k).eq.'F')cf=cf+1
		if(signstab(i,k).eq.'G')cg=cg+1
	  end do
	  nf=MAX(ca,cb,cc,cd,ce,cf,cg)
	  w=wcoef(nf,nseq)	 
	  write(8,801)npoz,nf,w				
	  wnp(npoz)=w
	end do
  end do
  close(8)

  do l=1,npoz
	if(l.gt.licz.OR.l.le.(npoz-licz))then
	  suma=0	
	  do k=l-licz,l+licz
		suma=suma+wnp(k)
	  end do
	  wind(k)=suma
	endif
  end do
  !	write(*,*)nres(1),nres(2),nres(3),nres(4),nres(5),nres(6)
  do l=1,npoz
	do i=1,nseq
	  m=nre(i,l)
	  !		write(*,*)nseq,l,i,m
	  if(m.ne.0)wnp2(i,m)=wnp(l)
	  if(m.ne.0)wind2(i,m)=wind(l)
	end do
  end do
  do i=1,nseq
	open(8,file=inpalign(1:(kp2-1))//'.wc.'//sfx(i),status='unknown',err=100)    
    do k=1,nres(i)
	  write(8,803)k,wind2(i,k)
	end do
	close(8)
  end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
elseif(algn_typ.eq.2)then
!...analyses the AA alignments
!.cut of the the gaps an numberizes motifs

!inicjalizacja licznikow		
	i=0
	nl=0
	nr=0
	ns=0

85	read(7,'(a)',end=86,err=99)line
	
if(line(1:1).ne.' ')then !ignoruje puste linie: linia z sekwencja w tym formacie zaczyna sie od opisu
nl=nl+1
lines(nl)=line !tabela "lines" ma tylko pelne linie 
i=i+1
!	kp=LEN(line)
!	if(i.eq.1)lseq=kp-1
else
	if(i.gt.0)then
	ns=i
	i=0
	endif
endif	
	
goto 85	

86 write(16,*)'alignment file read in'

OK: do kl=1,nl,ns
!write(*,*)'pocz',kl,nl,ns
	do j=29,78
		do i=1,ns
		if(lines(kl+i-1)(j:j).eq.' ')then
		EXIT OK
		endif
		end do
	if(lines(kl)(j:j).ne.' ')then
	nr=nr+1
		do i=1,ns
	signstab(i,nr)=lines(kl+i-1)(j:j)
		end do
!	write(*,*)ns,kl,signstab(1,nr),signstab(2,nr)
	endif
	end do
end do OK

!a=LEN(lines)
!b=LEN(lines)
!write(*,*)'dddddddd',a,b

!looks for gaps and cut them off
nmot=1
motl=0
do i=1,ns
numrs(i)=1
beg(i)=1
end do
nul=.FALSE.
full=.FALSE.
j=0

!write(*,*)'heeehheee',nr
SEK: do j=1,nr
!write(*,*)'plum',j
	do i=1,ns
	if(signstab(i,j).eq.'-')then
!	write(*,*)'chuj',j,i
	nul=.TRUE.
	endif
	end do
!write(*,*)j,nul	
	if(nul)then
		do i=1,ns
		
		if(signstab(i,j).ne.'-')then	
!jezeli nie jestes w bloku ale dana sekwencja nie ma akurat przerwy
!write(*,*)'sss',numrs(i),full
!			if(nmot.eq.1)then	!dla 1 motywu
!			if(.NOT.full)beg(i)=numrs(i)+1
!jezeli jeszcze sie nie zaczal 			
!			endif
		full=.FALSE.
		numrs(i)=numrs(i)+1
!	write(*,*)'beg',nmot,beg(i)
		endif
		end do
	endif

	if(nul)then
!	write(*,*)'chuj',j,motl
	do i=1,ns
		if(motl.eq.0)then
		nul=.FALSE.
		CYCLE SEK
		endif
	!gap found stop the counter
!		if(nmot.eq.1)then
!		beg(i)=
!		endif
	seq_motif(nmot,i,1)=beg(i)
	seq_motif(nmot,i,2)=kon(i)
!	write(*,*)'eee',j,i,nmot,beg(i),kon(i)
!	beg(i)=numrs(i)+1
	
	end do
	motl=0
	nul=.FALSE.
	nmot=nmot+1
	CYCLE SEK
	endif	

	do i=1,ns
	if(signstab(i,j).ne.'-')then
	if(.NOT.full)beg(i)=numrs(i)	
	kon(i)=numrs(i)
!	write(*,*)'test',j,i,beg(i),kon(i)
	numrs(i)=numrs(i)+1
		if(i.eq.1)motl=motl+1
!write(*,*)'test',j,i,numrs(i)
		endif
	end do
	full=.TRUE.
end do SEK
nmot=nmot-1
!write(*,*)'ssssss',nmot
!jezeli sekwencja konczy sie blokiem i nie ma odcinka z myslnikami
if(full)then
!write(*,*)'chuj',nul
	nmot=nmot+1
	do i=1,ns
	seq_motif(nmot,i,1)=beg(i)
	seq_motif(nmot,i,2)=kon(i)
!write(*,*)i,beg(i),kon(i)
	end do
endif

!if(numrs(1).lt.numrs(2))then
!maxprot=numrs(1)

!do i=1,ns
!if(numrs(i).gt.maxprot)then
!maxprot=numrs(i)
!endif
!end do

!write(*,*)'mot',seq_motif(1,1,1),seq_motif(1,1,2),seq_motif(1,2,1),seq_motif(1,2,2)

ignmtf=0
do i=1,ns
 do j=1,nmot
  if(bcutoff(1).gt.seq_motif(j,1,2).or.bcutoff(i).gt.seq_motif(j,i,2))then
  write(*,*)'WARNING: trim whole ',j,' motif'
  ignmtf=ignmtf+1
  endif
 end do
end do

!przesuwa numeracje o obciete z poczatku kawalki
!motyw 1: jezeli zaczyna od jedynki to zostaw
!write(*,*)'dddddddd',bcutoff(1),bcutoff(2),ecutoff(1)
if(seq_motif(1,1,1).eq.1)then
!jezeli bcutoff>0 skroceniu ulega pierwszy motyw: skroc tak samo 2 sekwencje
seq_motif2(1,1,1)=1
seq_motif2(1,1,2)=seq_motif(1,1,2)-bcutoff(1)
do i=1,ns
seq_motif2(1,i,1)=seq_motif(1,i,1)-bcutoff(i)+bcutoff(1)
seq_motif2(1,i,2)=seq_motif(1,i,2)-bcutoff(i)
end do
!else
!seq_motif2(1,1,1)=seq_motif(1,1,1)-bcutoff(1)
!seq_motif2(1,1,2)=seq_motif(1,1,2)-bcutoff(1)
endif

do i=1,ns
if(seq_motif(1,i,1).eq.1)then
!jezeli bcutoff>0 skroceniu ulega drugi motyw: skroc tak samo 1 sekwencje
seq_motif2(1,i,1)=1
seq_motif2(1,i,2)=seq_motif(1,i,2)-bcutoff(i)
seq_motif2(1,1,1)=seq_motif(1,1,1)-bcutoff(1)+bcutoff(i)
seq_motif2(1,1,2)=seq_motif(1,1,2)-bcutoff(1)
!else
!seq_motif2(1,numf,1)=seq_motif(1,numf,1)-bcutoff(numf)
!seq_motif2(1,numf,2)=seq_motif(1,numf,2)-bcutoff(numf)
endif
end do

!jezeli numer AA poczatku motywu jest mniejszy od 1 toz zrob =1 dla danej sekwencji a dla drugiej obetnij tyle ile trzeba :)
	if(seq_motif2(1,1,1).lt.1)then
	do i=1,ns
	seq_motif2(1,i,1)=seq_motif2(1,i,1)+ABS(seq_motif2(1,1,1))+1
	end do
	seq_motif2(1,1,1)=1
	endif
	do i=1,ns
	if(seq_motif2(1,i,1).lt.1)then
	seq_motif2(1,1,1)=seq_motif2(1,1,1)+ABS(seq_motif2(1,i,1))+1
	seq_motif2(1,i,1)=1
	endif
	end do

!write(*,*)1,seq_motif(1,1,1),seq_motif(1,1,2),seq_motif2(1,1,1),seq_motif2(1,1,2)
!write(*,*)1,seq_motif(1,numf,1),seq_motif(1,numf,2),seq_motif2(1,numf,1),seq_motif2(1,numf,2)

!dla pozostalych motywow po prostu odejmij

do j=2,nmot
seq_motif2(j,1,1)=seq_motif(j,1,1)-bcutoff(1)
seq_motif2(j,1,2)=seq_motif(j,1,2)-bcutoff(1)
do i=1,ns
seq_motif2(j,i,1)=seq_motif(j,i,1)-bcutoff(i)
seq_motif2(j,i,2)=seq_motif(j,i,2)-bcutoff(i)
end do
	if(seq_motif2(j,1,1).lt.1)then
!jezeli numer AA poczatku motywu jest mniejszy od 1 toz zrob =1 dla danej sekwencji a dla drugiej obetnij tyle ile trzeba :)
	do i=1,ns
	seq_motif2(j,i,1)=seq_motif2(j,i,1)+ABS(seq_motif2(j,1,1))+1
	end do
	seq_motif2(j,1,1)=1
	endif
	do i=1,ns
	if(seq_motif2(j,i,1).lt.1)then
	seq_motif2(j,1,1)=seq_motif2(j,1,1)+ABS(seq_motif2(j,i,1))+1
	seq_motif2(j,i,1)=1
	endif
	end do
!write(*,*)j,seq_motif(j,1,1),seq_motif(j,1,2),seq_motif2(j,1,1),seq_motif2(j,1,2)
!write(*,*)j,seq_motif(j,i,1),seq_motif(j,i,2),seq_motif2(j,i,1),seq_motif2(j,i,2)
end do

!motyw oststni: przytnij jesli potrzeba
do i=1,ns
seq_motif2(nmot,1,2)=seq_motif2(nmot,1,2)-ecutoff(1)-ecutoff(i)
seq_motif2(nmot,i,2)=seq_motif2(nmot,i,2)-ecutoff(i)-ecutoff(1)
end do

else
write(*,850)'algn_typ'
write(16,850)'algn_typ'
endif

goto 101

801	format(i5,1x,i2,1x,f6.3)
802	format(i85)
803	format(i5,1x,f6.3)
850     format('ERROR: Wrong value of parameter: ',a8,/&
       'restart with proper value set')
99	write(16,*)'ERROR reading alignment file:',inpalign
100   write(16,*)'ERROR opening unit: 7 or 8'

101	end
	
	
