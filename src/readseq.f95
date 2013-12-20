SUBROUTINE readseq(resnam_str,resid_str,nres_str,resnam_aa,resid_aa,nres_aa,nswitch,nasseq,aasekw,assekw)
USE tablice	
IMPLICIT none

character::line*90,comment*30,resnam_str(20,1000)*1,resnam_aa(20,1000)*3,aasekw*50,assekw*50
INTEGER::resid_str(20,1000),resid_aa(20,1000),nres_aa(20),nres_str(20),nswitch,nasseq,naaseq,n,m
		
!..reads in sequences in one letter format

!===========================================================
!..FASTA (Pearson) format

!-----------------------------------------------------------------
!...read in aminoacid sequence	
if(nswitch.eq.2.OR.nswitch.eq.4)then		
	open(7,file=aasekw,status='old',err=100)
	n=0
	naaseq=0
5	read(7,'(a)',end=42,err=100)line
	 if(line(1:1).eq.'>')then
	 comment=line(1:30)
	 naaseq=naaseq+1
	  if(naaseq.gt.1)then
	  nres_aa(naaseq-1)=n
	  write(16,*)'AA seq #: ',naaseq-1, ' lenght: ',n
	  endif
	 n=0
	 goto 5
	 elseif(line(1:1).eq.' ')then
	 goto 5
	 endif
	
	do 30 m=1,90
	 do 20 k=1,20
	 if(line(m:m).eq.AA(k))then
	 n=n+1
	 resnam_aa(naaseq,n)=aminok(k)
	 resid_aa(naaseq,n)=n
!	 write(*,*)naaseq,n,resnam_aa(naaseq,n)
	 goto 30
	 endif
20	 continue
30	continue
	goto 5
	close(7)
endif
42	nres_aa(naaseq)=n	
write(16,*)'AA seq #: ',naaseq, ' lenght: ',nres_aa(naaseq)
write(16,*)'Number of AA sequences: ',naaseq	
!---------------------------------------------------------------
!...reads in structural sequence
if(nswitch.eq.3.OR.nswitch.eq.4)then		
  open(7,file=assekw,status='old',err=100)
  n=0
  nasseq=0
50	read(7,'(a)',end=82,err=100)line
	if(line(1:1).eq.'>')then
	comment=line(1:30)
	nasseq=nasseq+1
	 if(nasseq.gt.1)then
	 nres_str(nasseq-1)=n
	 write(16,*)'SA seq #: ',nasseq-1, ' lenght: ',n
	 endif
	n=0
	goto 50
	elseif(line(1:1).eq.' ')then
	goto 50
	endif
	
	do 70 m=1,90
	 do 60 k=1,7
	 if(line(m:m).eq.AS(k))then
	 n=n+1
	 resnam_str(nasseq,n)=AS(k)
	 resid_str(nasseq,n)=n
!	 write(*,*)nasseq,n,resnam_str(nasseq,n)
	 goto 70
	 endif
60	 continue
70	continue	
	goto 50
close(7)
endif
82	nres_str(nasseq)=n
write(16,*)'SA seq #: ',nasseq, ' lenght: ',nres_str(nasseq)
write(16,*)'Number of SA sequences: ',nasseq	

if(naaseq.ne.nasseq)then
  write(16,*)'WARNING: AA and SA sequence numbers not equal!'
  write(16,*)'Be sure that thera are same numbers of AA and SA seuqences in the input files' 
endif

goto 101	
	
	
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
	
100	write(16,*)'ERR opening unit:',16
101	end
		