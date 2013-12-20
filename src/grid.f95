SUBROUTINE grid(crd,nres,nseq,kwezel,space,bscenter,ucenter,ngcat,kcenter,grdtab,atomy,bocz,resid,&
grdcon,seed,agrid,gpk,gpn,gridtype,ligand,nlig)
IMPLICIT none
!**********************  creates the grid  *****************************
  
INTEGER::bocz(15000),gpk(-30:30,-30:30,-30:30,2),agrid(-30:30,-30:30,-30:30,2),nres,nseq,ngcat,kcenter,&
gpn(-30:30,-30:30,-30:30,2),kwezel(3),cgrd,resid(15000),gridtype,nlig,kwezx,kwezy,kwezz,c1count,&
 c2count,c3count,c4count,na,na2,mg,mgw,i,k,l,m,n,lres,lat,nwezx,nwezy,nwezz,nr,nr2
REAL(8)::crd(100000,3,2),atomy(15000,30,4),seed(3),space(3),bscenter(3),ucenter(3),grdtab(15000,30),&
grdcon(15000,30),ligand(100,3),shels(4,100),px,py,pz,cx,cy,cz,coord1,coord2,coord3,coord4,odst,sumx,&
sumy,sumz,xh,yh,zh,odl
     
     
!###################################################################################
	
	if(gridtype.eq.1.OR.gridtype.eq.3)then
!..grid type 1: binding site of the second molecule is oriented to 
!..overlay the ligands. grid is created for two molecule in the same
!..origin based on the first grid
!
!..grid type 3: siatki sa porownane bez nakladania ligandow - trzeba uwazac 
!..aby centrum dla obu siatek bylo porownywalne np. centrum geometryczne liganda	
	
	
!...defines grid center
!...possible grid centers: center of the coordinate space, 
!...geometrical center of a molecule, geometric center of 
!...the binding site, user defined coordinates, defined atom 
	
	if(gridtype.eq.1.AND.nseq.ne.1)goto 15
	if(kcenter.eq.0.)then
!...origin of the system	
	seed(1)=0.0
	seed(2)=0.0
	seed(3)=0.0
	
	elseif(kcenter.eq.1)then
!...defined atom coordinates	
	 do 3 i=1,3
	 seed(i)=INT(crd(ngcat,i,nseq)/space(i))*space(i)        
3	 continue

	elseif(kcenter.eq.2)then
!...geometric center	
!...defines geometric center of a molecule/chain	
	sumx=0
	sumy=0
	sumz=0
	lat=0
	do 7 k=1,nres
		do 6 n=1,bocz(k)	
	       sumx=sumx+atomy(k,n,2)
	       sumy=sumy+atomy(k,n,3)
	       sumz=sumz+atomy(k,n,4)
	       lat=lat+1
6		continue
7	continue
	seed(1)=sumx/lat
	seed(2)=sumy/lat
	seed(3)=sumz/lat

	elseif(kcenter.eq.3)then
!...center of the binding site	
	seed(1)=bscenter(1)
	seed(2)=bscenter(2)
	seed(3)=bscenter(3)
		
	elseif(kcenter.eq.4)then
!...user defined coordinates	
	seed(1)=ucenter(1)
	seed(2)=ucenter(2)
	seed(3)=ucenter(3)
	
	else
	write(16,850)'kcenter'
	stop
	endif
	
	do 13 l=-kwezel(1),kwezel(1)
	 do 12 m=-kwezel(2),kwezel(2)
	  do 11 n=-kwezel(3),kwezel(3)
	  agrid(l,m,n,1)=0
	  gpk(l,m,n,1)=0
	  gpn(l,m,n,1)=0
!	  write(16,*)l,m,n,agrid(l,m,n,1),agrid(l,m,n,2)
11	  continue
12	 continue
13	continue
	

!	write(*,*)seed(1),seed(2),seed(3),bscenter(1),bscenter(2),bscenter(3),nres	
!	write(*,*)kwezel(1)-1,kwezel(2)-1,kwezel(3)-1,nseq   	
	
15	if(nseq.eq.2)then
	do 18 l=-kwezel(1),kwezel(1)
	 do 17 m=-kwezel(2),kwezel(2)
	  do 16 n=-kwezel(3),kwezel(3)
	  agrid(l,m,n,2)=0
!	  write(16,*)l,m,n,agrid(l,m,n,1),agrid(l,m,n,2)
16	  continue
17	 continue
18	continue
	endif
	
      do 30 k=1,nres
       do 29 n=1,bocz(k)
!...defining grid points    
	cx=atomy(k,n,2)-seed(1)       
	cy=atomy(k,n,3)-seed(2)
	cz=atomy(k,n,4)-seed(3)
        nwezx=AINT(cx/space(1))        
        nwezy=AINT(cy/space(2))
        nwezz=AINT(cz/space(3))
	if(ABS(nwezx).gt.kwezel(1)-1.OR.ABS(nwezy).gt.kwezel(2)-1.OR.ABS(nwezz).gt.kwezel(3)-1)then
	goto 29
	endif
	if(nwezx.lt.0)nwezx=nwezx-1
	if(nwezy.lt.0)nwezy=nwezy-1
	if(nwezz.lt.0)nwezz=nwezz-1
        px=(nwezx+0.5)*space(1)
        py=(nwezy+0.5)*space(2)
        pz=(nwezz+0.5)*space(3)
!	write(16,*)k,n,cx,nwezx,px,cy,nwezy,py,cz,nwezz,pz
        if(cx.le.px)then
         if(cy.le.py)then
          if(cz.le.pz)then
!..point # 1          
		gpk(nwezx,nwezy,nwezz,nseq)=k
		gpn(nwezx,nwezy,nwezz,nseq)=n
		agrid(nwezx,nwezy,nwezz,nseq)=1
	  else
!..point # 5 	  
		gpk(nwezx,nwezy,nwezz+1,nseq)=k
		gpn(nwezx,nwezy,nwezz+1,nseq)=n
		agrid(nwezx,nwezy,nwezz+1,nseq)=1
	  endif	
	 else
	  if(cz.le.pz)then
!..point # 3          
	      	gpk(nwezx,nwezy+1,nwezz,nseq)=k
	      	gpn(nwezx,nwezy+1,nwezz,nseq)=n
	        agrid(nwezx,nwezy+1,nwezz,nseq)=1
	  else     	
!..point # 7 
	  	gpk(nwezx,nwezy+1,nwezz+1,nseq)=k
	  	gpn(nwezx,nwezy+1,nwezz+1,nseq)=n
	  	agrid(nwezx,nwezy+1,nwezz+1,nseq)=1
	  endif
	 endif 
	else
	 if(cy.le.py)then
	  if(cz.le.pz)then
!..point # 2          
	        gpk(nwezx+1,nwezy,nwezz,nseq)=k
	        gpn(nwezx+1,nwezy,nwezz,nseq)=n
	        agrid(nwezx+1,nwezy,nwezz,nseq)=1
	  else
!..point # 6
		gpk(nwezx+1,nwezy,nwezz+1,nseq)=k
		gpn(nwezx+1,nwezy,nwezz+1,nseq)=n
		agrid(nwezx+1,nwezy,nwezz+1,nseq)=1
	  endif	
	 else
	  if(cz.le.pz)then
!..point # 4          
		gpk(nwezx+1,nwezy+1,nwezz,nseq)=k
		gpn(nwezx+1,nwezy+1,nwezz,nseq)=n
		agrid(nwezx+1,nwezy+1,nwezz,nseq)=1
	  else     	
!..point # 8 
		gpk(nwezx+1,nwezy+1,nwezz+1,nseq)=k
		gpn(nwezx+1,nwezy+1,nwezz+1,nseq)=n
		agrid(nwezx+1,nwezy+1,nwezz+1,nseq)=1
	  endif
	 endif 
	endif
!...washing up old table
	grdtab(k,n)=0	
29	 continue
30	continue
	
	if(nseq.eq.2)then
!	write(*,*)'aaaaa',kwezel(1),kwezel(2),kwezel(3)	
!...writes consensus grid points values        
	write(4,811)"x node","y node","z node","AA #","AA id","atom #","node value"      	
	do 45 l=-kwezel(1),kwezel(1)
	 do 42 m=-kwezel(2),kwezel(2)
	  do 40 n=-kwezel(3),kwezel(3) 
	  mgw=agrid(l,m,n,1)
	  mg=agrid(l,m,n,2)
	  nr2=gpk(l,m,n,2)
	  na2=gpn(l,m,n,2)
	  nr=gpk(l,m,n,1)
	  na=gpn(l,m,n,1)
!	  write(16,*)l,m,n,mgw,mg,nr2,na2
	  if(mgw.eq.0.AND.mg.eq.0)cgrd=0
	  if(mgw.eq.1.AND.mg.eq.0)cgrd=1
	  if(mgw.eq.0.AND.mg.eq.1)cgrd=2
	  if(mgw.eq.1.AND.mg.eq.1)cgrd=3
	  grdtab(nr2,na2)=cgrd
!	  write(*,*)'dsfsfds',nr,na,grdtab(nr,na)
	   if(cgrd.eq.3)then
	   grdcon(nr,na)=grdcon(nr,na)+1
	   endif
	  lres=resid(nr)	
	  write(4,810)l,m,n,nr,lres,na,cgrd
40	  continue
42	 continue
45	continue
	endif
!	write(*,*)'dsfsfds',grdtab(188,2)
	
	goto 100
	
!###################################################################################   	
	elseif(gridtype.eq.2)then
	
!..grid type 2: siatki niezalezne, zliczane kontakty w czterech strefach 
!..koordynacyjnych dla kazdego liganda	
	
	
!...defines grid centers for each molecule separately
!...possible grid centers: center of the coordinate space, 
!...geometrical center of a molecule, geometric center of 
!...the binding site, user defined coordinates, defined atom 
	
	
	if(kcenter.eq.0.)then
!...origin of the system	
	seed(1)=0.0
	seed(2)=0.0
	seed(3)=0.0
	
	elseif(kcenter.eq.1)then
!...defined atom coordinates	
	 do 53 i=1,3
	 seed(i)=INT(crd(ngcat,i,nseq)/space(i))*space(i)        
53	 continue

	elseif(kcenter.eq.2)then
!...geometric center	
!...defines geometric center of a molecule/chain	
	sumx=0
	sumy=0
	sumz=0
	lat=0
	do 57 k=1,nres
		do 56 n=1,bocz(k)	
	       sumx=sumx+atomy(k,n,2)
	       sumy=sumy+atomy(k,n,3)
	       sumz=sumz+atomy(k,n,4)
	       lat=lat+1
56		continue
57	continue
	seed(1)=sumx/lat
	seed(2)=sumy/lat
	seed(3)=sumz/lat

	elseif(kcenter.eq.3)then
!...center of the binding site	
	seed(1)=bscenter(1)
	seed(2)=bscenter(2)
	seed(3)=bscenter(3)
		
	elseif(kcenter.eq.4)then
!...user defined coordinates	
	seed(1)=ucenter(1)
	seed(2)=ucenter(2)
	seed(3)=ucenter(3)
	
	else
	write(16,850)'kcenter'
	stop
	endif
	
	if(nseq.eq.1)then
	do 65 l=-kwezel(1),kwezel(1)
	 do 64 m=-kwezel(2),kwezel(2)
	  do 63 n=-kwezel(3),kwezel(3)
	  agrid(l,m,n,1)=0
	  gpk(l,m,n,1)=0
	  gpn(l,m,n,1)=0
!	  write(16,*)l,m,n,agrid(l,m,n,1),agrid(l,m,n,2)
63	  continue
64	 continue
65	continue
	endif

!	write(*,*)seed(1),seed(2),seed(3),bscenter(1),bscenter(2),bscenter(3),nres	
!	write(*,*)kwezel(1)-1,kwezel(2)-1,kwezel(3)-1,nseq   	
	
	if(nseq.eq.2)then
	do 68 l=-kwezel(1),kwezel(1)
	 do 67 m=-kwezel(2),kwezel(2)
	  do 66 n=-kwezel(3),kwezel(3)
	  agrid(l,m,n,2)=0
!	  write(16,*)l,m,n,agrid(l,m,n,1),agrid(l,m,n,2)
66	  continue
67	 continue
68	continue
	endif
	
      do 70 k=1,nres
       do 69 n=1,bocz(k)
!...defining grid points    
	cx=atomy(k,n,2)-seed(1)       
	cy=atomy(k,n,3)-seed(2)
	cz=atomy(k,n,4)-seed(3)
        nwezx=AINT(cx/space(1))        
        nwezy=AINT(cy/space(2))
        nwezz=AINT(cz/space(3))
	if(ABS(nwezx).gt.kwezel(1)-1.OR.ABS(nwezy).gt.kwezel(2)-1.OR.ABS(nwezz).gt.kwezel(3)-1)then
	goto 69
	endif
	if(nwezx.lt.0)nwezx=nwezx-1
	if(nwezy.lt.0)nwezy=nwezy-1
	if(nwezz.lt.0)nwezz=nwezz-1
        px=(nwezx+0.5)*space(1)
        py=(nwezy+0.5)*space(2)
        pz=(nwezz+0.5)*space(3)
!	write(16,*)k,n,cx,nwezx,px,cy,nwezy,py,cz,nwezz,pz
        if(cx.le.px)then
         if(cy.le.py)then
          if(cz.le.pz)then
!..point # 1          
		gpk(nwezx,nwezy,nwezz,nseq)=k
		gpn(nwezx,nwezy,nwezz,nseq)=n
		agrid(nwezx,nwezy,nwezz,nseq)=1
	  else
!..point # 5 	  
		gpk(nwezx,nwezy,nwezz+1,nseq)=k
		gpn(nwezx,nwezy,nwezz+1,nseq)=n
		agrid(nwezx,nwezy,nwezz+1,nseq)=1
	  endif	
	 else
	  if(cz.le.pz)then
!..point # 3          
	      	gpk(nwezx,nwezy+1,nwezz,nseq)=k
	      	gpn(nwezx,nwezy+1,nwezz,nseq)=n
	        agrid(nwezx,nwezy+1,nwezz,nseq)=1
	  else     	
!..point # 7 
	  	gpk(nwezx,nwezy+1,nwezz+1,nseq)=k
	  	gpn(nwezx,nwezy+1,nwezz+1,nseq)=n
	  	agrid(nwezx,nwezy+1,nwezz+1,nseq)=1
	  endif
	 endif 
	else
	 if(cy.le.py)then
	  if(cz.le.pz)then
!..point # 2          
	        gpk(nwezx+1,nwezy,nwezz,nseq)=k
	        gpn(nwezx+1,nwezy,nwezz,nseq)=n
	        agrid(nwezx+1,nwezy,nwezz,nseq)=1
	  else
!..point # 6
		gpk(nwezx+1,nwezy,nwezz+1,nseq)=k
		gpn(nwezx+1,nwezy,nwezz+1,nseq)=n
		agrid(nwezx+1,nwezy,nwezz+1,nseq)=1
	  endif	
	 else
	  if(cz.le.pz)then
!..point # 4          
		gpk(nwezx+1,nwezy+1,nwezz,nseq)=k
		gpn(nwezx+1,nwezy+1,nwezz,nseq)=n
		agrid(nwezx+1,nwezy+1,nwezz,nseq)=1
	  else     	
!..point # 8 
		gpk(nwezx+1,nwezy+1,nwezz+1,nseq)=k
		gpn(nwezx+1,nwezy+1,nwezz+1,nseq)=n
		agrid(nwezx+1,nwezy+1,nwezz+1,nseq)=1
	  endif
	 endif 
	endif
!...washing up old table
	grdtab(k,n)=0	
69	 continue
70	continue

!..orienting ligand


	do 95 i=1,nlig
	xh=ligand(i,1)-seed(1)
	yh=ligand(i,2)-seed(2)
	zh=ligand(i,3)-seed(3)
	nwezx=AINT(xh/space(1))        
        nwezy=AINT(yh/space(2))
        nwezz=AINT(zh/space(3))
	if(ABS(nwezx).gt.kwezel(1)-1.OR.ABS(nwezy).gt.kwezel(2)-1.OR.ABS(nwezz).gt.kwezel(3)-1)then
	write(*,*)'WARNING: ligand outside the receptor grid.' 
     	write(*,*)'Increase the grid size'
	endif
			
!	write(*,*)'aaaaa',kwezel(1),kwezel(2),kwezel(3)	
!...writes consensus grid points values        
!	write(2,811)"x node","y node","z node","AA #","AA id","atom #",
!     &	"node value"

     	kwezx=AINT(12/space(1))
	kwezy=AINT(12/space(2))
	kwezz=AINT(12/space(3))
	
	odst=MIN(space(1),space(2),space(3))
	
	coord1=2.5/odst
	coord2=5/odst
	coord3=7.5/odst
	coord4=10/odst
	c1count=0
	c2count=0
	c3count=0
	c4count=0	
	
	do 85 l=nwezx-kwezx,nwezx+kwezx
	 do 82 m=nwezy-kwezy,nwezy+kwezy
	  do 80 n=nwezz-kwezz,nwezz+kwezz
	  
	  odl=SQRT((space(l)**2)+(space(m)**2)+(space(n)**2))
	  if(odl.lt.coord1)then
	  c1count=c1count+1
	  elseif(odl.lt.coord2)then
	  c2count=c2count+1
	  elseif(odl.lt.coord3)then
	  c3count=c3count+1
	  elseif(odl.lt.coord4)then
	  c4count=c4count+1
	  endif
80	  continue
82	 continue
85	continue

	shels(1,i)=c1count
	shels(2,i)=c2count
	shels(3,i)=c3count
	shels(4,i)=c4count
95	continue

	endif
!	write(*,*)'dsfsfds',grdtab(188,2)
	
	goto 100
        
	
		
810     format(3(8x,I3),3(6x,I5),10x,I1) 
811	format(3(5x,a6),7x,a4,5x,a5,6x,a6,1x,a10)

850   format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
	write(16,*)'grid computing...OK'
100      end
