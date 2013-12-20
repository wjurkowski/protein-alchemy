SUBROUTINE volume(nazwa,elliptic_vol,bocz,atomy,am,volref,nres,ellatomy)
IMPLICIT none
    
REAL(8)::v,ellatomy(15000,30,4),xcent(15000),ycent(15000),zcent(15000),xmcent(15000),ymcent(15000),zmcent(15000),&
atomy(15000,30,4),am(15000,30),at(15000,30,4),a,b,c,buf,sumx,sumy,sumz,RM,vch,xmin,xmax,ymin,ymax,zmin,zmax
character*15 nazwa
INTEGER::bocz(15000),elliptic_vol,volref,nres,nat,k,l,m,n

if(elliptic_vol.eq.1)then
  do k=1,nres
	do n=1,bocz(k)
	  at(k,n,2)=ellatomy(k,n,2)
	  at(k,n,3)=ellatomy(k,n,3)
	  at(k,n,4)=ellatomy(k,n,4)
	end do
  end do
elseif(elliptic_vol.eq.0)then
  do k=1,nres
	do n=1,bocz(k)
	  at(k,n,2)=atomy(k,n,2)
	  at(k,n,3)=atomy(k,n,3)
	  at(k,n,4)=atomy(k,n,4)
	end do
  end do
else
  write(16,850)'elliptic_vol'
endif

if(volref.eq.0)then
	xmin=at(1,2,2)
	ymin=at(1,2,3)
	zmin=at(1,2,4)
	xmax=at(1,2,2)
	ymax=at(1,2,3)
	zmax=at(1,2,4)
	do k=1,nres
!	write(*,*)k,xmax,xmin
      if(at(k,2,2).lt.xmin)xmin=at(k,2,2)
      if(at(k,2,2).gt.xmax)xmax=at(k,2,2)
      if(at(k,2,3).lt.ymin)ymin=at(k,2,3)
      if(at(k,2,3).gt.ymax)ymax=at(k,2,3)
      if(at(k,2,4).lt.zmin)zmin=at(k,2,4)
      if(at(k,2,4).gt.zmax)zmax=at(k,2,4)
	end do  

elseif(volref.eq.1)then
	xmin=at(1,1,2)
	ymin=at(1,1,3)
	zmin=at(1,1,4)
	xmax=at(1,1,2)
	ymax=at(1,1,3)
	zmax=at(1,1,4)
	do k=1,nres
	  do n=1,bocz(k)
		if(at(k,n,2).lt.xmin)xmin=at(k,n,2)
		if(at(k,n,2).gt.xmax)xmax=at(k,n,2)
		if(at(k,n,3).lt.ymin)ymin=at(k,n,3)
		if(at(k,n,3).gt.ymax)ymax=at(k,n,3)
		if(at(k,n,4).lt.zmin)zmin=at(k,n,4)
		if(at(k,n,4).gt.zmax)zmax=at(k,n,4)
	  end do
	end do
elseif(volref.eq.2)then
	!...geometric center of a side chain	
	do k=1,nres
	  sumx=0
	  sumy=0
	  sumz=0
	  do n=5,bocz(k)	
        sumx=sumx+at(k,n,2)
        sumy=sumy+at(k,n,3)
        sumz=sumz+at(k,n,4)
	  end do
	  xcent(k)=sumx/(n-4)
	  ycent(k)=sumy/(n-4)
	  zcent(k)=sumz/(n-4)
	end do
    xmin=xcent(1)
    xmax=xcent(1)
    ymin=ycent(1)
    ymax=ycent(1)
    zmin=zcent(1)
    zmax=zcent(1)
    do k=1,nres
      if(xcent(k).lt.xmin)xmin=xcent(k)
      if(xcent(k).gt.xmax)xmax=xcent(k)
      if(ycent(k).lt.ymin)ymin=ycent(k)
      if(ycent(k).gt.ymax)ymax=ycent(k)
      if(zcent(k).lt.zmin)zmin=zcent(k)
      if(zcent(k).gt.zmax)zmax=zcent(k)
	end do
	!write(*,*)'xmin', xmin,'xmax',xmax
elseif(volref.eq.3)then
	!...mass center of a side chain	
	do k=1,nres
	  sumx=0
	  sumy=0
	  sumz=0
	  RM=0.0
	  do n=5,bocz(k)	
        sumx=sumx+at(k,n,2)*am(k,n)
        sumy=sumy+at(k,n,3)*am(k,n)
        sumz=sumz+at(k,n,4)*am(k,n)
		RM=RM+am(k,n)
	  end do
	  xmcent(k)=sumx/RM
	  ymcent(k)=sumy/RM
	  zmcent(k)=sumz/RM
	end do
    xmin=xmcent(1)
    xmax=xmcent(1)
    ymin=ymcent(1)
    ymax=ymcent(1)
    zmin=zmcent(1)
    zmax=zmcent(1)
    do k=1,nres
      if(xmcent(k).lt.xmin)xmin=xmcent(k)
      if(xmcent(k).gt.xmax)xmax=xmcent(k)
      if(ymcent(k).lt.ymin)ymin=ymcent(k)
      if(ymcent(k).gt.ymax)ymax=ymcent(k)
      if(zmcent(k).lt.zmin)zmin=zmcent(k)
      if(zmcent(k).gt.zmax)zmax=zmcent(k)	
	end do
elseif(volref.eq.4)then
	!...geometric center of a residue	
	do k=1,nres
	  sumx=0
	  sumy=0
	  sumz=0
	  do n=1,bocz(k)	
        sumx=sumx+at(k,n,2)
        sumy=sumy+at(k,n,3)
        sumz=sumz+at(k,n,4)
	  end do
	  xcent(k)=sumx/n
	  ycent(k)=sumy/n
	  zcent(k)=sumz/n
	end do
    xmin=xcent(1)
    xmax=xcent(1)
    ymin=ycent(1)
    ymax=ycent(1)
    zmin=zcent(1)
    zmax=zcent(1)
    do k=1,nres
      if(xcent(k).lt.xmin)xmin=xcent(k)
      if(xcent(k).gt.xmax)xmax=xcent(k)
      if(ycent(k).lt.ymin)ymin=ycent(k)
      if(ycent(k).gt.ymax)ymax=ycent(k)
      if(zcent(k).lt.zmin)zmin=zcent(k)
      if(zcent(k).gt.zmax)zmax=zcent(k)
	end do
elseif(volref.eq.5)then
	!...center of mass of a residue	
	do k=1,nres
	  sumx=0
	  sumy=0
	  sumz=0
	  RM=0.0
	  do n=1,bocz(k)	
        sumx=sumx+at(k,n,2)*am(k,n)
        sumy=sumy+at(k,n,3)*am(k,n)
        sumz=sumz+at(k,n,4)*am(k,n)
		RM=RM+am(k,n)
	  end do
	  xmcent(k)=sumx/RM
	  ymcent(k)=sumy/RM
	  zmcent(k)=sumz/RM
	end do
    xmin=xmcent(1)
    xmax=xmcent(1)
    ymin=ymcent(1)
    ymax=ymcent(1)
    zmin=zmcent(1)
    zmax=zmcent(1)
    do k=1,nres
      if(xmcent(k).lt.xmin)xmin=xmcent(k)
      if(xmcent(k).gt.xmax)xmax=xmcent(k)
      if(ymcent(k).lt.ymin)ymin=ymcent(k)
      if(ymcent(k).gt.ymax)ymax=ymcent(k)
      if(zmcent(k).lt.zmin)zmin=zmcent(k)
      if(zmcent(k).gt.zmax)zmax=zmcent(k)	
	end do
else
	write(16,850)'volref'
    stop
endif
  
a=ABS(xmax-xmin)
b=ABS(ymax-ymin)
c=ABS(zmax-zmin)
v=a*b*c
if(a.gt.b)then
	buf=b
	b=a
	a=buf
endif
write(16,*)'writing output file...'
write(*,'(a10,1x,i4,1x,3f8.3,1x,f14.3)')nazwa,nres,a,b,c,v 
!write(*,*)nazwa,nres,a,b,c,v
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
      
end
