SUBROUTINE dihsort(fipsi,grmesh,ngrcv,nres,aadihsave)
USE tablice
IMPLICIT none

REAL(8)::fipsi(15000,3),grmesh
INTEGER::ngrcv(-180:180,-180:180,20),kltot(20),ngrfi,ngrpsi,ifi,ipsi,llp,m,nres,aadihsave
character*15 outfile,line*80,text*10,text2,cfile*12,key*1,wsad,ans*1,nazwa
	
      do 15 i=1,20
       do 10 ifi=-180,180,INT(grmesh)
        do 5 ipsi=-180,180,INT(grmesh)
       ngrcv(ifi,ipsi,i)=0
5       continue
10     continue
15    continue

	do 16 i=1,20	
	kltot(i)=0
16	continue

        llp=llp+1
	
if(aadihsave.eq.1)then	
	do 18 i=1,20
        outfile=(aminok(i)//'.s')
      open(20+i,file=outfile,status='unknown')
18	continue
        
26     do 28 m=1,nres
	do 27 i=1,20
	if(fipsi(m,3).eq.i)then
        write(20+i,810)fipsi(m,1),fipsi(m,2)
	endif
27	continue
28     continue   
810    format(2f8.3) 

	do 29 i=1,20
	endfile 20+i      
       close(20+i,status='keep')
29	continue
elseif(aadihsave.ne.0)then
write(*,850)'aadihsave'
stop
endif

!..counting dihedrals values sorted out with defined step


	do 55 m=1,nres
         do 50 i=1,20
	if(fipsi(m,3).eq.i)then
       ngrfi=(INT(fipsi(m,1)/grmesh))*grmesh
       ngrpsi=(INT(fipsi(m,2)/grmesh))*grmesh
	if(ngrfi.lt.0.)ngrfi=ngrfi-INT(grmesh)
	if(ngrpsi.lt.0.)ngrpsi=ngrpsi-INT(grmesh)
        ngrcv(ngrfi,ngrpsi,i)=ngrcv(ngrfi,ngrpsi,i)+1
	endif
50       continue
55	continue

850   format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
       
       end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
