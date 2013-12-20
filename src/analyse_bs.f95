SUBROUTINE analyse_bs(atomy,coord,calfc,nres,cutoff_bs,resid,residue_ref,nseq,am,bocz,asseq,&
het_present,ligand,nlig,resnam,nmol,resmap,bscenter,motdl,namcor,bsalign,fnamcor,aaseq,savefit,ncach,nch,&
molbm)
USE tablice
IMPLICIT none

REAL(8)::atomy(15000,30,4),coord(100000,3,2),ligand(100,3),xcent(15000),ycent(15000),zcent(15000),&
xmcent(15000),ymcent(15000),zmcent(15000),am(15000,30),calfc(15000,3,2),motif(7,4,3,99,1000),&
bscenter(3),Xm(3,28),Ym(3,28),Xpt(3,28),rmsd,d,sumx,sumy,sumz,xsum,ysum,zsum,RM,x,y,z,xh,yh,zh
INTEGER::bocz(15000),resid(15000),cmap(100,15000),resmap(15000),molbm(1000),motres(1000,99),cutoff_bs,&
ncach(0:50),nres,cnter,residue_ref,nseq,nlig,nmol,motdl,bsalign,savefit,nch,ltab,lpos,kpos,j,kk,l,m,n,&
lat,kres
character*4 asseq(15000)*1,motyw1,motyw2,motyw3,motyw4,namcor*50,aaseq(15000)*1,resnam(15000)*3,&
fnamcor*50,name*50
logical het_present,contact
		
	kpos=index(namcor,' ')
	lpos=index(fnamcor,' ')
	
	if(.NOT.het_present)then
	write(*,*)'hetero atoms must be present!'
	goto 200
	endif

!...aa codes translation
        do 3 n=1,nres
         do 2 l=1,20
        if(resnam(n).eq.aminok(l))then
        aaseq(n)=AA(l)
        endif
2        continue
3       continue

!...looks for residues close to ligand within a defined cutoff
!...residue can be tought as an effective atom  
	do 80 i=1,nlig
	
	xh=ligand(i,1)
	yh=ligand(i,2)
	zh=ligand(i,3)

	if(residue_ref.eq.0)then
      do 5 k=1,nres
      x=calfc(k,1,nseq)
      y=calfc(k,2,nseq)
      z=calfc(k,3,nseq)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	if(d.le.cutoff_bs)then
	cmap(i,k)=1
	else
	cmap(i,k)=0
	endif
5     continue

   	elseif(residue_ref.eq.1)then
      do 15 k=1,nres
      	contact=.FALSE.
      	do 10 n=1,bocz(k)
	x=atomy(k,n,2)
        y=atomy(k,n,3)
        z=atomy(k,n,4)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	if(d.le.cutoff_bs)then
	cmap(i,k)=1
	goto 15
	else  
	cmap(i,k)=0
	endif
10      continue	
15    continue     
	
   	elseif(residue_ref.eq.2)then
!...geometric center of a side chain	
	do 17 k=1,nres
	sumx=0
	sumy=0
	sumz=0
		do 16 n=5,bocz(k)	
        sumx=sumx+atomy(k,n,2)
        sumy=sumy+atomy(k,n,3)
        sumz=sumz+atomy(k,n,4)
16		continue
	xcent(k)=sumx/(n-4)
	ycent(k)=sumy/(n-4)
	zcent(k)=sumz/(n-4)
17	continue
      do 25 k=1,nres
      x=xcent(k)
      y=ycent(k)
      z=zcent(k)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	if(d.le.cutoff_bs)then
	cmap(i,k)=1
	else
	cmap(i,k)=0
	endif
25    continue     

   	elseif(residue_ref.eq.3)then
!...mass center of a side chain	
	do 27 k=1,nres
	sumx=0
	sumy=0
	sumz=0
	RM=0.0
		do 26 n=5,bocz(k)	
        sumx=sumx+atomy(k,n,2)*am(k,n)
        sumy=sumy+atomy(k,n,3)*am(k,n)
        sumz=sumz+atomy(k,n,4)*am(k,n)
	RM=RM+am(k,n)
26		continue
	xmcent(k)=sumx/RM
	ymcent(k)=sumy/RM
	zmcent(k)=sumz/RM
27	continue
      do 35 k=1,nres
      x=xmcent(k)
      y=ymcent(k)
      z=zmcent(k)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	if(d.le.cutoff_bs)then
	cmap(i,k)=1
	else
	cmap(i,k)=0
	endif
35	continue
	
   	elseif(residue_ref.eq.4)then
!...geometric center of a residue	
	do 37 k=1,nres
	sumx=0
	sumy=0
	sumz=0
		do 36 n=1,bocz(k)	
        sumx=sumx+atomy(k,n,2)
        sumy=sumy+atomy(k,n,3)
        sumz=sumz+atomy(k,n,4)
36		continue
	xcent(k)=sumx/n
	ycent(k)=sumy/n
	zcent(k)=sumz/n
37	continue
      do 45 k=1,nres
      x=xcent(k)
      y=ycent(k)
      z=zcent(k)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	if(d.le.cutoff_bs)then
	cmap(i,k)=1
	else
	cmap(i,k)=0
	endif
45    continue     
 
  	elseif(residue_ref.eq.5)then
!...center of mass of a residue	
	do 47 k=1,nres
	sumx=0
	sumy=0
	sumz=0
	RM=0.0
		do 46 n=1,bocz(k)	
        sumx=sumx+atomy(k,n,2)*am(k,n)
        sumy=sumy+atomy(k,n,3)*am(k,n)
        sumz=sumz+atomy(k,n,4)*am(k,n)
	RM=RM+am(k,n)
46		continue
	xmcent(k)=sumx/RM
	ymcent(k)=sumy/RM
	zmcent(k)=sumz/RM
47	continue

	do 55 k=1,nres
      	x=xmcent(k)
      	y=ymcent(k)
      	z=zmcent(k)
	d=SQRT((x-xh)**2+(y-yh)**2+(z-zh)**2)
	 if(d.le.cutoff_bs)then
	 cmap(i,k)=1
	 else
	 cmap(i,k)=0
	 endif
55	continue	

	else
      	write(16,850)'residue_ref'
      	stop
      	endif
		
80	continue
	
	
	do 83 k=1,nres
	do 82 i=1,nlig
	if(cmap(i,k).eq.1)then
	resmap(k)=1
	goto 83
	endif	
82	continue
	resmap(k)=0
83	continue	
	
	
	cnter=0
	xsum=0.0
	ysum=0.0
	zsum=0.0
	lat=0
	
!...saving binding site structure codes	
	write(12,809)"res #","res id","pos4","pos3","pos2","pos1"
	write(13,809)"res #","res id","pos4","pos3","pos2","pos1"
        
	do 90 k=1,nres
	if(resmap(k).eq.1)then
	
	do 85 n=1,bocz(k)
	xsum=xsum+atomy(k,n,2)
	ysum=ysum+atomy(k,n,3)
	zsum=zsum+atomy(k,n,4)
	lat=lat+1
85	continue	
	
	 if(k.ge.4)then
	  if(k.le.nres-3)then
	 motyw1=asseq(k-3)//asseq(k-2)//asseq(k-1)//asseq(k)
	 motyw2=asseq(k-2)//asseq(k-1)//asseq(k)//asseq(k+1)
	 motyw3=asseq(k-1)//asseq(k)//asseq(k+1)//asseq(k+2)
	 motyw4=asseq(k)//asseq(k+1)//asseq(k+2)//asseq(k+3)
        write(12,804)k,resid(k),motyw1,motyw2,motyw3,motyw4
	 motyw1=aaseq(k-3)//aaseq(k-2)//aaseq(k-1)//aaseq(k)
	 motyw2=aaseq(k-2)//aaseq(k-1)//aaseq(k)//aaseq(k+1)
	 motyw3=aaseq(k-1)//aaseq(k)//aaseq(k+1)//aaseq(k+2)
	 motyw4=aaseq(k)//aaseq(k+1)//aaseq(k+2)//aaseq(k+3)
        write(13,804)k,resid(k),motyw1,motyw2,motyw3,motyw4
	
	cnter=cnter+1
	      do 88 l=1,motdl	
		do 87 n=1,4
		motif(l,n,1,cnter,nmol)=atomy(k-4+l,n,2)
		motif(l,n,2,cnter,nmol)=atomy(k-4+l,n,3)
		motif(l,n,3,cnter,nmol)=atomy(k-4+l,n,4)
87		continue
88	      continue	
	motres(nmol,cnter)=k

	  elseif(k.eq.nres-2)then
	 motyw1=asseq(k-3)//asseq(k-2)//asseq(k-1)//asseq(k)
	 motyw2=asseq(k-2)//asseq(k-1)//asseq(k)//asseq(k+1)
	 motyw3=asseq(k-1)//asseq(k)//asseq(k+1)//asseq(k+2)
        write(12,803)k,resid(k),motyw1,motyw2,motyw3
	 motyw1=aaseq(k-3)//aaseq(k-2)//aaseq(k-1)//aaseq(k)
	 motyw2=aaseq(k-2)//aaseq(k-1)//aaseq(k)//aaseq(k+1)
	 motyw3=aaseq(k-1)//aaseq(k)//aaseq(k+1)//aaseq(k+2)
        write(13,803)k,resid(k),motyw1,motyw2,motyw3
	  elseif(k.eq.nres-1)then
	 motyw1=asseq(k-3)//asseq(k-2)//asseq(k-1)//asseq(k)
	 motyw2=asseq(k-2)//asseq(k-1)//asseq(k)//asseq(k+1)
        write(12,802)k,resid(k),motyw1,motyw2
	 motyw1=aaseq(k-3)//aaseq(k-2)//aaseq(k-1)//aaseq(k)
	 motyw2=aaseq(k-2)//aaseq(k-1)//aaseq(k)//aaseq(k+1)
        write(13,802)k,resid(k),motyw1,motyw2
	  elseif(k.eq.nres)then
	 motyw1=asseq(k-3)//asseq(k-2)//asseq(k-1)//asseq(k)
        write(12,801)k,resid(k),motyw1
	 motyw1=aaseq(k-3)//aaseq(k-2)//aaseq(k-1)//aaseq(k)
        write(13,801)k,resid(k),motyw1
	  endif	
	 elseif(k.eq.3)then
	 motyw2=asseq(k-2)//asseq(k-1)//asseq(k)//asseq(k+1)
	 motyw3=asseq(k-1)//asseq(k)//asseq(k+1)//asseq(k+2)
	 motyw4=asseq(k)//asseq(k+1)//asseq(k+2)//asseq(k+3)
        write(12,805)k,resid(k),motyw2,motyw3,motyw4
	 motyw2=aaseq(k-2)//aaseq(k-1)//aaseq(k)//aaseq(k+1)
	 motyw3=aaseq(k-1)//aaseq(k)//aaseq(k+1)//aaseq(k+2)
	 motyw4=aaseq(k)//aaseq(k+1)//aaseq(k+2)//aaseq(k+3)
        write(13,805)k,resid(k),motyw2,motyw3,motyw4
	 elseif(k.eq.2)then
	 motyw3=asseq(k-1)//asseq(k)//asseq(k+1)//asseq(k+2)
	 motyw4=asseq(k)//asseq(k+1)//asseq(k+2)//asseq(k+3)
        write(12,806)k,resid(k),motyw3,motyw4
	 motyw3=aaseq(k-1)//aaseq(k)//aaseq(k+1)//aaseq(k+2)
	 motyw4=aaseq(k)//aaseq(k+1)//aaseq(k+2)//aaseq(k+3)
        write(13,806)k,resid(k),motyw3,motyw4
	 elseif(k.eq.1)then
	 motyw4=asseq(k)//asseq(k+1)//asseq(k+2)//asseq(k+3)
        write(12,807)k,resid(k),motyw4
	 motyw4=aaseq(k)//aaseq(k+1)//aaseq(k+2)//aaseq(k+3)
        write(13,807)k,resid(k),motyw4
	 endif
	
	endif
90      continue

	molbm(nmol)=cnter
	bscenter(1)=xsum/lat
	bscenter(2)=ysum/lat
	bscenter(3)=zsum/lat
	
	write(4,812)"res id","contact"
	do 92 k=1,nres
	write(4,811)resid(k),resmap(k)
92	continue

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	
!....superimpose structural motifes from a set of chains
if(bsalign.eq.1)then
  if(nmol.ge.2)then
	 ltab=motdl*4
!..aligns binding-site motifes of the same ligand
!...alignes pairwaise first motifs with the rest
	name=fnamcor(1:(lpos-1))//"."//sfx(nmol)//namcor(1:(kpos-1))//'.bsrmsd'
	open(20+nmol,file=name,status='unknown',access='APPEND',action='readwrite')
	write(20+nmol,814)'pair #','mol # 1','mol # ','RMSD'
     
     kk=0
	 do 140 j=1,molbm(1)
	  name=namcor(1:(kpos-1))//"."//sfx(j)//'.bsan.pdb'
	  if(savefit.eq.1)then
	  open(3,file=name,status='unknown')
	  elseif(savefit.ne.0)then
	  write(16,850)'savefit'
	  stop
	  endif
	 
	 m=0
 	      do 104 l=1,motdl
                do 103 n=1,4
                 m=m+1
                 Ym(1,m)=motif(l,n,1,j,1)
                 Ym(2,m)=motif(l,n,2,j,1)
                 Ym(3,m)=motif(l,n,3,j,1)	
103               continue
104             continue

	do 130 i=1,molbm(nmol)
		m=0
              do 112 l=1,motdl
                do 111 n=1,4
		m=m+1
                Xm(1,m)=motif(l,n,1,i,nmol)
                Xm(2,m)=motif(l,n,2,i,nmol)
                Xm(3,m)=motif(l,n,3,i,nmol)
111             continue
112           continue
	
	call fitting(Xm,Ym,ltab,rmsd)

	if(savefit.eq.1)then
!...saves new motifes coordinates
	m=0
	kres=motres(nmol,i)
	do 125 l=1,motdl
	 m=m+1
	 write(3,813)'ATOM',m,'N',resnam(kres-4+l),i,resid(kres-4+l),Xpt(1,m),Xpt(2,m),Xpt(3,m)
	 m=m+1
	 write(3,813)'ATOM',m,'CA',resnam(kres-4+l),i,resid(kres-4+l),Xpt(1,m),Xpt(2,m),Xpt(3,m)
 	m=m+1
	 write(3,813)'ATOM',m,'C',resnam(kres-4+l),i,resid(kres-4+l),Xpt(1,m),Xpt(2,m),Xpt(3,m)
	 m=m+1
	 write(3,813)'ATOM',m,'O',resnam(kres-4+l),i,resid(kres-4+l),Xpt(1,m),Xpt(2,m),Xpt(3,m)	
125	continue
	 write(3,'(a)')'TER'
	elseif(savefit.ne.0)then
	write(16,850)'savefit'
	stop
	endif
	
!...saves rmsd values for given motif
	kk=kk+1
	write(20+nmol,815)kk,j,i,rmsd
	
130	  continue
	if(savefit.eq.1)close(3)			
140       continue
  close(20+nmol,status='keep')
  endif
elseif(bsalign.ne.0)then
  write(16,850)'bsalign'
  stop
endif
				
	
200	write(16,*)'...done'
	

801     format(i7,1x,i7,1x,1x,a4)
802     format(i7,1x,i7,1x,1x,a4,1x,a4)
803     format(i7,1x,i7,1x,1x,a4,1x,a4,1x,a4)
804     format(i7,1x,i7,1x,1x,a4,1x,a4,1x,a4,1x,a4)
805     format(i7,1x,i7,1x,6x,a4,1x,a4,1x,a4)
806     format(i7,1x,i7,1x,11x,a4,1x,a4)
807     format(i7,1x,i7,1x,16x,a4)
809	format(a7,1x,a7,1x,1x,a4,1x,a4,1x,a4,1x,a4)
811	format(i5,1x,i1)
812	format(a7,1x,a7)
813	format(a4,1x,i6,2x,a2,2x,a3,i2,1x,i3,4x,3f8.3,24x)
814	format(a10,1x,a10,1x,a10,1x,a10)
815	format(7x,i3,8x,i3,8x,i3,1x,f10.6)
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')

end



