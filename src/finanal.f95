SUBROUTINE finanal(ngrcv,grmesh,dihstat,grdcon,gridok,fnamcor,flin1,flin2,fatom,fbocz,fnres)
USE tablice
IMPLICIT none

INTEGER::ngrcv(-180:180,-180:180,20),kltot(20),kl(20),fnres,fbocz(15000),ifi,ipsi,kpos,n,totsum,gridok,&
dihstat
real(8)::pfipsi(-180:180,-180:180,20),fatom(15000,30,4),grdcon(15000,30),grmesh
character flin1(15000,30)*30,flin2(15000,30)*24,cfile*12,out*2,fnamcor*50
	
kpos=index(fnamcor,' ')

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx		
if(dihstat.eq.1)then
  write(16,*)'performs fi psi distribution analysis...'
  !.....performs statistical analysis of fi psi distribution
  !.....writes fi,psi values count
  do i=1,20
	kl(i)=0
    do ifi=-180,180,INT(grmesh)
	  do ipsi=-180,180,INT(grmesh)
		!	write(*,816)ifi,ipsi,ngrcv(ifi,ipsi,i)
		if(ngrcv(ifi,ipsi,i).ne.0)then
        kl(i)=kl(i)+ngrcv(ifi,ipsi,i)
		endif
	  end do
	end do
	kltot(i)=kltot(i)+kl(i)
  end do
!.................calculates probability
!...total sum
  totsum=0
  do i=1,20
	  totsum=totsum+kltot(i)
	  do ifi=-180,180,INT(grmesh)
		do ipsi=-180,180,INT(grmesh)
		  pfipsi(ifi,ipsi,i)=ngrcv(ifi,ipsi,i)/kltot(i)
		end do
	  end do
  end do
  write(16,*)'...distribution analysis completed'
elseif(dihstat.ne.0)then
  write(*,850)'dihstat'
  stop
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	
if(gridok.eq.1)then
  write(16,*)'...writes structure with consensus grid'	
  open(3,file=(fnamcor(1:(kpos-1))//'.consgr.pdb'),status='unknown',err=90)
  do k=1,fnres
	do n=1,fbocz(k)
	  write(3,802)flin1(k,n),fatom(k,n,2),fatom(k,n,3),fatom(k,n,4),flin2(k,n)(1:6),grdcon(k,n),flin2(k,n)(13:24)
	end do
  end do
  write(3,801)'TER'
  close(3)
elseif(gridok.ne.0)then
  write(*,850)'gridok'
  stop
endif	

goto 100
	
90 	write(16,*)'ERR opening unit 3'
stop		
		
801   format(a3)
802   format(a30,3f8.3,a6,f6.2,a12)	
816    format(I4,2X,I4,2X,I8)
818     format(' '/&
               'total:'/&
               a3,' = ',I8)
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')    
 
100	end
