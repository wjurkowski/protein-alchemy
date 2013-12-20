SUBROUTINE crd2pdb(crd,nseq,hetco,nonhet,ellatomy,bocz,linia1,linia2,nres,ellpdb,nch,ncach,namcor,motifpdb,resmap,&
atomy,grdtab,gridok,outpdb,minres,tfrepl,fatom,estrtype,aasekw,resnam_aa,resnam_str,nres_str)
USE tablice
IMPLICIT none
!.. writes output to the pdb format file

character*80 line,text1*30,text2*24,linia1(15000,30)*30,linia2(15000,30)*24,namcor*50,nazwa*16,restype*3,aasekw*50,&
resnam_aa(20,1000)*3,resnam_str(20,1000)*1
real(8)::crd(100000,3,2),hetco(10000,3,2),ellatomy(15000,30,4),atomy(15000,30,4),grdtab(15000,30),tfrepl(15000,30),&
fatom(15000,30,4)
INTEGER::bocz(15000),ncach(0:50),resmap(15000),nseq,nonhet,nch,nres,kpos,l,n,motifpdb,gridok,outpdb,minres,ellpdb,estrtype,&
lpos,nres_str(20)
kpos=index(namcor,' ')

!..xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if(outpdb.eq.1)then
  do l=1,nch
    if((ncach(l)-ncach(l-1)).lt.minres)cycle
	open(3,file=namcor(1:(kpos-1))//"."//sfx(l)//'.o.pdb',status='unknown',err=90)
	write(3,'(a)')'REMARK created with Siatami8.5'	
	do k=ncach(l-1)+1,ncach(l)
	  do n=1,bocz(k)
		write(3,800)linia1(k,n),atomy(k,n,2),atomy(k,n,3),atomy(k,n,4),tfrepl(k,n),linia2(k,n)(15:24)
	  end do
	end do	
	write(3,801)'TER'
	close(3)
  end do
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	   
if(ellpdb.eq.1)then
  if(estrtype.eq.1)then 
	write(16,*)'...writing elliptic structure based on read in native structure:'
	!.....writes elliptic structure based on native structure
	do l=1,nch
	  if((ncach(l)-ncach(l-1)).lt.minres)cycle
	  open(3,file=namcor(1:(kpos-1))//"."//sfx(l)//'.e.pdb',status='unknown',err=90)
	  do k=ncach(l-1)+1,ncach(l)
		do n=1,bocz(k)
		  write(3,800)linia1(k,n),ellatomy(k,n,2),ellatomy(k,n,3),ellatomy(k,n,4),linia2(k,n)
		end do
	  end do  
	  write(3,801)'TER'
	  close(3)
	end do
  elseif(estrtype.eq.2)then
	lpos=index(aasekw,' ')
	do i=1,1
	  open(3,file=aasekw(1:(lpos-1))//"."//sfx(i)//'.e.pdb',status='unknown',err=90)
	!....writes elliptic structure pdb based on read in structural sequence	
	!..begining tetraalanine
	l=0
	do k=1,4
	  write(3,804)'ATOM',l,'N','ALA',i,k,ellatomy(k,1,2),ellatomy(k,1,3),ellatomy(k,1,4)
      l=l+1
	  write(3,804)'ATOM',l,'CA','ALA',i,k,ellatomy(k,2,2),ellatomy(k,2,3),ellatomy(k,2,4)
      l=l+1
	  write(3,804)'ATOM',l,'C','ALA',i,k,ellatomy(k,3,2),ellatomy(k,3,3),ellatomy(k,3,4)
      l=l+1
	  write(3,804)'ATOM',l,'O','ALA',i,k,ellatomy(k,4,2),ellatomy(k,4,3),ellatomy(k,4,4)
	end do
	do k=5,nres_str(i)-4
	  l=l+1	
	  write(3,803)'ATOM',l,'N',resnam_aa(i,k-4),i,k,ellatomy(k,1,2),ellatomy(k,1,3),ellatomy(k,1,4),resnam_str(i,k-4)
      l=l+1
	  write(3,803)'ATOM',l,'CA',resnam_aa(i,k-4),i,k,ellatomy(k,2,2),ellatomy(k,2,3),ellatomy(k,2,4),resnam_str(i,k-4)
      l=l+1
	  write(3,803)'ATOM',l,'C',resnam_aa(i,k-4),i,k,ellatomy(k,3,2),ellatomy(k,3,3),ellatomy(k,3,4),resnam_str(i,k-4)
      l=l+1
	  write(3,803)'ATOM',l,'O',resnam_aa(i,k-4),i,k,ellatomy(k,4,2),ellatomy(k,4,3),ellatomy(k,4,4),resnam_str(i,k-4)
	end do	
	do k=nres_str(i)-3,nres_str(i)
	  write(3,804)'ATOM',l,'N','ALA',i,k,ellatomy(k,1,2),ellatomy(k,1,3),ellatomy(k,1,4)
      l=l+1
	  write(3,804)'ATOM',l,'CA','ALA',i,k,ellatomy(k,2,2),ellatomy(k,2,3),ellatomy(k,2,4)
      l=l+1
	  write(3,804)'ATOM',l,'C','ALA',i,k,ellatomy(k,3,2),ellatomy(k,3,3),ellatomy(k,3,4)
      l=l+1
	  write(3,804)'ATOM',l,'O','ALA',i,k,ellatomy(k,4,2),ellatomy(k,4,3),ellatomy(k,4,4)
	end do
	write(3,'(a)')'TER'

	close(3)
	end do

	write(16,*)'...done'
  endif
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if(motifpdb.eq.1)then
  write(16,*)'...writes binding-site motifes structure'			
  open(3,file=namcor(1:(kpos-1))//".bsall.pdb",status='unknown',err=91)
  do k=1,nres
	if(resmap(k).eq.1)then
	  if(k.ge.4)then
		if(k.le.nres-3)then
		  do l=1,7
			do n=1,4
			  write(3,800)linia1(k-4+l,n),atomy(k-4+l,n,2),atomy(k-4+l,n,3),atomy(k-4+l,n,4),linia2(k-4+l,n)
			end do
		  end do
		  write(3,801)'TER'
		endif
	  endif 
	endif	
  end do
  close(3)
endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	
if(gridok.eq.1)then
	if(nseq.eq.2)then
	  write(16,*)'...writes structure of grided protein'
	  open(3,file=(namcor(1:(kpos-1))//'.gr.pdb'),status='unknown',err=90)
	  do k=1,nres
		do n=1,bocz(k)
		  write(3,802)linia1(k,n),atomy(k,n,2),atomy(k,n,3),atomy(k,n,4),linia2(k,n)(1:6),grdtab(k,n),linia2(k,n)(13:24)
		end do
	  end do
	  write(3,801)'TER'
	  close(3)
	endif
endif

goto 100
800   format(a30,3f8.3,5x,f8.3,1x,a10)
802   format(a30,3f8.3,a6,f6.2,a12)
801   format(a3)
803   format(a30,3f8.3,a6,a24)
804   format(a4,1x,i6,2x,a2,2x,a3,i2,1x,i3,4x,3f8.3,3x)
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')

90 	  write(16,*)'ERR opening unit 4'
91 	  write(16,*)'ERR opening unit 45'
100   end
