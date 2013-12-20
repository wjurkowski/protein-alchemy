SUBROUTINE ramadist(fipsi,ellfipsi,nch,ncach,minres,namcor,resid,resnam,eceppok,emesh,kwsize2,tfrepl,bocz,nres,meanr)
USE tablice
IMPLICIT none

real(8)::ellfipsi(15000,2),fipsi(15000,3),enertab(20,360,360),dist1(nres),dist2(nres),dist3(nres),dist4(nres),dist7(nres),&
dist8(nres),dist9(nres),dist10(nres),dist11(nres),dist12(nres),dist13(nres),dist14(nres),tfrepl(15000,30),fi,psi,efi,epsi,&
Eell,En,emesh,suma(12),dcart,dmanh,dczeb,dcanb
integer ncach(0:50),resid(15000),bocz(15000),enpsi,enepsi,enfi,enefi,nch,minres,kwsize2,kpos,licz,lami,l,m,n,nres,eceppok,meanr
character*1 resnam(15000)*3,namcor*50,tab*1

!-----------------------------------------------------------------------	
!......reading ecepp energy maps
if(eceppok.eq.1)then
	write(16,*)'reads ecepp based energy maps...'
	call eceppanal(emesh,enertab)
	write(16,*)'...maps read in'
elseif(eceppok.ne.0)then
      	write(*,850)'eceppok'
       	stop
endif
tab=char(9)
kpos=index(namcor,' ')
licz=INT(kwsize2/2)
!.....calculates native-elliptic fi, psi distances
write(16,*)'calculates native-elliptic fi, psi distances...'

do m=1,nch
  if((ncach(m)-ncach(m-1)).lt.minres)cycle
  open(18,file=namcor(1:(kpos-1))//"."//sfx(m)//".rdist",status='unknown')
  if(meanr.eq.0)then
	if(eceppok.eq.0)write(18,809)"AA_#",tab,"AA_id",tab,"AA",tab,"dcart",tab,"dmanh",tab,"dczeb",tab,"dcanb",tab,"dfi",tab,"dpsi"
!	if(eceppok.eq.0)write(18,809)"AA_#","AA_id","AA","dcart","dmanh","dczeb","dcanb","dfi","dpsi"
	if(eceppok.eq.1)write(18,810)"AA_#",tab,"AA_id",tab,"AA",tab,"dcart",tab,"dmanh",tab,"dczeb",tab,"dcanb",tab,"dfi",tab,&
	"dpsi",tab,"dE",tab,"dE/dcart"
  end if 
	do k=ncach(m-1)+1,ncach(m)-1
	  do l=1,20
		if(aminok(l).eq.resnam(k))lami=l
	  end do
	  efi=ellfipsi(k,1)
      epsi=ellfipsi(k,2)
      fi=fipsi(k,1)
      psi=fipsi(k,2)

	  dist1(k)=dcart(fi,efi,psi,epsi)!..euklides..
	  dist2(k)=dmanh(fi,efi,psi,epsi)!..manhattan
	  dist3(k)=dczeb(fi,efi,psi,epsi)!..czebujew
	  dist4(k)=dcanb(fi,efi,psi,epsi)!..canberr
	  dist7(k)=efi-fi!fi difference
	  dist8(k)=epsi-psi!psi difference
	  if(eceppok.eq.1)then
		enfi=AINT((fi+180)/emesh)+1
		enefi=AINT((efi+180)/emesh)+1
		enpsi=AINT((psi+180)/emesh)+1
		enepsi=AINT((epsi+180)/emesh)+1
		En=enertab(lami,enpsi,enfi)
		Eell=enertab(lami,enepsi,enefi)
		dist9(k)=En-Eell
		dist10(k)=dist9(k)/dist1(k)
	  endif
	  if(eceppok.eq.0)write(18,812)k,tab,resid(k),tab,resnam(k),tab,dist1(k),tab,dist2(k),tab,dist3(k),tab,dist4(k),tab,dist7(k),&
	  tab,dist8(k)
	  if(eceppok.eq.1)write(18,813)k,tab,resid(k),tab,resnam(k),tab,dist1(k),tab,dist2(k),tab,dist3(k),tab,dist4(k),tab,dist7(k),&
	  dist8(k),tab,dist9(k),tab,dist10(k)
	end do

	if(meanr.eq.1)then
	  if(eceppok.eq.0)write(18,815)"AA_#",tab,"AA_id",tab,"    AA",tab,"dcart_w",tab,"dmanh_w",tab,"dczeb_w",tab,"dcanb_w",tab,&
	  "dfi",tab,"dpsi"
	  if(eceppok.eq.1)write(18,811)"AA_#",tab,"AA_id",tab,"AA",tab,"dcart_w",tab,"dmanh_w",tab,"dczeb_w",tab,"dcanb_w",tab,&
	  "dfi_w",tab,"dpsi_w",tab,"dE_w",tab,"dE/dcart_w",tab,"|dE|_w"
	  do k=ncach(m-1)+1,ncach(m)-1
		if(k.gt.licz.OR.k.le.(ncach(m)-1-licz))then
		  do i=1,12
			suma(i)=0
		  end do
		  do l=k-licz,k+licz
			suma(1)=suma(1)+dist1(l)
			suma(2)=suma(2)+dist2(l)
			suma(3)=suma(3)+dist3(l)
			suma(4)=suma(4)+dist4(l)
			suma(7)=suma(7)+dist7(l)
			suma(8)=suma(8)+dist8(l)
			if(eceppok.eq.1)then
			  suma(9)=suma(9)+dist9(l)
			  suma(10)=suma(10)+dist10(l)
			  suma(11)=suma(11)+ABS(dist9(l))
			  suma(12)=suma(12)+ABS(dist10(l))
			endif
		  end do
		  dist1(k)=suma(1)/kwsize2
		  dist2(k)=suma(2)/kwsize2
		  dist3(k)=suma(3)/kwsize2
		  dist4(k)=suma(4)/kwsize2
		  dist7(k)=suma(7)/kwsize2
		  dist8(k)=suma(8)/kwsize2	
		  if(eceppok.eq.0)write(18,812)k,tab,resid(k),tab,resnam(k),tab,dist1(k),tab,dist2(k),tab,dist3(k),tab,dist4(k),tab,&
		  dist7(k),tab,dist8(k)
		  if(eceppok.eq.1)then
			dist9(k)=suma(9)/kwsize2
			dist10(k)=suma(10)/kwsize2
			dist11(k)=suma(11)/kwsize2
			write(18,814)k,tab,resid(k),tab,resnam(k),tab,dist1(k),tab,dist2(k),tab,dist3(k),tab,dist4(k),tab,dist7(k),tab,&
			dist8(k),tab,dist9(k),tab,dist10(k),tab,dist11(k)
		endif 
	  endif
	  !..puts calculated value as beta factor
	  do n=1,bocz(k)
		tfrepl(k,n)=dist1(k)
	  end do
	end do
	endif
  close(18)
end do

809		format(9(A8,A1))
810		format(11(A10,A1))
811		format(12(A10,A1))
812		format(2(I8,A1),A8,A1,6(F8.3,A1))
813		format(2(I8,A1),A8,A1,8(F8.3,A1))	
814		format(2(I8,A1),A8,A1,9(F8.3,A1))
815		format(9(A8,A1))
850		format('ERROR: Wrong value of parameter: ',a8,/&
       'restart with proper value set')
END
