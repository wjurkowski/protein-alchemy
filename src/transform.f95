SUBROUTINE transform(coord,calfc,hetco,nhet,lp,lpca,nseq,refsys,planecrd,nch,ncach,natch,atomy,bocz,nonhet)
IMPLICIT none
       
REAL(8)::cootr(100000,3,2),x(2),y(2),z(2),space(3),hetco(10000,3,2),calfc(15000,3,2),coord(100000,3,2),&
xt(100000),yt(100000),zt(100000),atomy(15000,30,4),planecrd(4,3,2),TRM(3,3),trans(3),A(4,4),B(4,3),dxy,&
dxz,xold,yold,zold,yx1,yx2,yy1,yy2,zx1,zx2,zy1,zy2,zz1,zz2,yx2po,yy2po,zx2po,zy2po,zx2po2,zz2po2,&
dist,cosfi,coste,sinte,sinfi,d,ymax,zmax
integer lpca(2),lp(2),grdco(15000,3),bocz(15000),natch(0:50),ncach(0:50),nhet,nseq,refsys,nch,nnch,LDA,LDB,NRHS,N,nonhet,mmax,&
ll,IPIV,INFO,i,j,k,l,m,kmax
EXTERNAL DGESV
      
!.......uklad wspolrzednych
3     if(refsys.eq.1)then
	goto 10 
      elseif(refsys.eq.2)then
      goto 40
      else
      write(16,850)'refsys'
      stop
      endif

    
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 10   N=4
      NRHS=3
      LDA=4
      LDB=4	
	
        do 12 i=1,4
       A(i,1)=planecrd(i,1,2)  
       A(i,2)=planecrd(i,2,2)  
       A(i,3)=planecrd(i,3,2)
       A(i,4)=1.
!       write(*,*)'A: ',A(i,1),A(i,2),A(i,3),A(i,4)
12    continue

      do 19 j=1,4 
        do 18 i=1,3       
       B(j,i)=planecrd(j,i,1)
18      continue
!	write(*,*)'B: ',B(j,1),B(j,2),B(j,3)
19    continue
              
      call DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
	
      do 22 j=1,3
            do 21 i=1,3
	TRM(i,j)=b(i,j)
21          continue
	trans(j)=b(4,j)
!	write(*,*)B(1,j),B(2,j),B(3,j),B(4,j)
!	write(*,*)TRM(1,j),TRM(2,j),TRM(3,j),trans(j)
22    continue


!..   transform structure using previous calculated transform matrix
	do 26 j=1,nch
      	do 25 i=natch(j-1)+1,natch(j)
      	xold=coord(i,1,2)
	yold=coord(i,2,2)
      	zold=coord(i,3,2)
    coord(i,1,nseq)=xold*TRM(1,1)+yold*TRM(2,1)+zold*TRM(3,1)+trans(1)
	coord(i,2,nseq)=xold*TRM(1,2)+yold*TRM(2,2)+zold*TRM(3,2)+trans(2)
	coord(i,3,nseq)=xold*TRM(1,3)+yold*TRM(2,3)+zold*TRM(3,3)+trans(3)
25    	continue
26    	continue	

	do 30 j=1,nch
      	do 29 k=ncach(j-1)+1,ncach(j)
	 do 28 n=1,bocz(k)
      	xold=atomy(k,n,2)
	yold=atomy(k,n,3)
      	zold=atomy(k,n,4)
	atomy(k,n,2)=xold*TRM(1,1)+yold*TRM(2,1)+zold*TRM(3,1)+trans(1)
	atomy(k,n,3)=xold*TRM(1,2)+yold*TRM(2,2)+zold*TRM(3,2)+trans(2)
      atomy(k,n,4)=xold*TRM(1,3)+yold*TRM(2,3)+zold*TRM(3,3)+trans(3)
28	 continue	
29    	continue
30    	continue	

	if(nonhet.eq.0)then
	do 35 k=1,nhet
	xold=hetco(k,1,2)
	yold=hetco(k,2,2)
	zold=hetco(k,3,2)
   	hetco(k,1,2)=xold*TRM(1,1)+yold*TRM(2,1)+zold*TRM(3,1)+trans(1)
	hetco(k,2,2)=xold*TRM(1,2)+yold*TRM(2,2)+zold*TRM(3,2)+trans(2)
	hetco(k,3,2)=xold*TRM(1,3)+yold*TRM(2,3)+zold*TRM(3,3)+trans(3)
!	write(*,*)hetco(k,1,2),hetco(k,2,2),hetco(k,3,2)
35	continue
	elseif(nonhet.ne.1)then
	write(16,850)'nonhet'
	endif

	goto 100	

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...orientacja osi-z zgodnie z najwieksza odlegloscia c-alfa w czasteczce    
40    ll=nch
      nnch=1
      

      do 80 l=1,nnch
     
	zmax=0.0
	ymax=0.0
      do 45 k=ncach(l-1)+1,ncach(ll)-1
                do 42 m=k+1,ncach(ll)
        x(1)=calfc(k,1,nseq)
        y(1)=calfc(k,2,nseq)
        z(1)=calfc(k,3,nseq)
        x(2)=calfc(m,1,nseq)
        y(2)=calfc(m,2,nseq)
        z(2)=calfc(m,3,nseq)
        dist=SQRT((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
        if(dist.gt.zmax)then
                zmax=dist
                kmax=k
                mmax=m
            if(z(1).lt.z(2))then
            zx1=x(1)
              zx2=x(2)
            zy1=y(1)
              zy2=y(2)
            zz1=z(1)
              zz2=z(2)
            else
            zx1=x(2)
              zx2=x(1)
              zy1=y(2)
              zy2=y(1)
              zz1=z(2)
              zz2=z(1)        
            endif
      endif
42              continue
45      continue
      
!...macierz obrotu
      write(16,*)'zmax ',zmax
      write(16,*)'zx2: ',zx2,' zy2: ',zy2,' zz2: ',zz2
      zx2=zx2-zx1
      zy2=zy2-zy1
      zz2=zz2-zz1
      write(16,*)'zx2: ',zx2,' zy2: ',zy2,' zz2: ',zz2
      
      dxy=SQRT(zx2**2+zy2**2)
      sinfi=zy2/dxy
      cosfi=zx2/dxy
      zx2po=zx2*cosfi+zy2*sinfi
      zy2po=zx2*(-sinfi)+zy2*cosfi
      if(ABS(zx2po).lt.9.0E-6)zx2po=0.0 
      if(ABS(zy2po).lt.9.0E-6)zy2po=0.0      
      write(16,*)'zx2: ',zx2po,' zy2: ',zy2po,' zz2: ',zz2     
      
      dxz=SQRT(zx2po**2+zz2**2)
      sinte=zx2po/dxz
      coste=zz2/dxz    
      zx2po2=zx2po*coste+zz2*(-sinte)
      zz2po2=zx2po*sinte+zz2*coste
      if(ABS(zx2po2).lt.9.0E-6)zx2po2=0.0 
      if(ABS(zz2po2).lt.9.0E-6)zz2po2=0.0      
      write(16,*)'zx2: ',zx2po2,' zy2: ',zy2po,' zz2: ',zz2po2  
      
!...obroty
        do 49 k=natch(l-1)+1,natch(ll)
        xt(k)=coord(k,1,nseq)-zx1
        yt(k)=coord(k,2,nseq)-zy1
        zt(k)=coord(k,3,nseq)-zz1
                  
      coord(k,1,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*coste+zt(k)*(-sinte)
      coord(k,2,nseq)=xt(k)*(-sinfi)+yt(k)*cosfi
      coord(k,3,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*sinte+zt(k)*coste
      if(ABS(coord(k,1,nseq)).lt.9.0E-6)coord(k,1,nseq)=0.0
      if(ABS(coord(k,2,nseq)).lt.9.0E-6)coord(k,2,nseq)=0.0
      if(ABS(coord(k,3,nseq)).lt.9.0E-6)coord(k,3,nseq)=0.0
49      continue

        do 51 k=ncach(l-1)+1,ncach(ll)
		do 50 n=1,bocz(k)
        xt(k)=atomy(k,n,2)-zx1
        yt(k)=atomy(k,n,3)-zy1
        zt(k)=atomy(k,n,4)-zz1
                  
      atomy(k,n,2)=(xt(k)*cosfi+yt(k)*sinfi)*coste+zt(k)*(-sinte)
      atomy(k,n,3)=xt(k)*(-sinfi)+yt(k)*cosfi
      atomy(k,n,4)=(xt(k)*cosfi+yt(k)*sinfi)*sinte+zt(k)*coste
      if(ABS(atomy(k,n,2)).lt.9.0E-6)atomy(k,n,2)=0.0
      if(ABS(atomy(k,n,3)).lt.9.0E-6)atomy(k,n,3)=0.0
      if(ABS(atomy(k,n,4)).lt.9.0E-6)atomy(k,n,4)=0.0
50		continue
51      continue

      do 55 k=ncach(l-1)+1,ncach(ll)
        xt(k)=calfc(k,1,nseq)-zx1
        yt(k)=calfc(k,2,nseq)-zy1
        zt(k)=calfc(k,3,nseq)-zz1
             
      calfc(k,1,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*coste+zt(k)*(-sinte)
      calfc(k,2,nseq)=xt(k)*(-sinfi)+yt(k)*cosfi
      calfc(k,3,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*sinte+zt(k)*coste
      if(ABS(calfc(k,1,nseq)).lt.9.0E-6)calfc(k,1,nseq)=0.0
      if(ABS(calfc(k,2,nseq)).lt.9.0E-6)calfc(k,2,nseq)=0.0
      if(ABS(calfc(k,3,nseq)).lt.9.0E-6)calfc(k,3,nseq)=0.0
55    continue

	if(nonhet.eq.0)then
      do 60 k=1,nhet
        xt(k)=hetco(k,1,nseq)-zx1
        yt(k)=hetco(k,2,nseq)-zy1
        zt(k)=hetco(k,3,nseq)-zz1
      
      hetco(k,1,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*coste+zt(k)*(-sinte)
      hetco(k,2,nseq)=xt(k)*(-sinfi)+yt(k)*cosfi
      hetco(k,3,nseq)=(xt(k)*cosfi+yt(k)*sinfi)*sinte+zt(k)*coste     
      if(ABS(hetco(k,1,nseq)).lt.(9.0E-6))hetco(k,1,nseq)=0.0
      if(ABS(hetco(k,2,nseq)).lt.(9.0E-6))hetco(k,2,nseq)=0.0
      if(ABS(hetco(k,3,nseq)).lt.(9.0E-6))hetco(k,3,nseq)=0.0
60    continue
	elseif(nonhet.ne.1)then
	write(16,850)'nonhet'
	endif

!...orientacja osi y wzgledem najwiekszej odleglosci c-alfa     
      do 68 k=ncach(l-1)+1,ncach(ll)-1
            do 67 m=k+1,ncach(ll)
        x(1)=calfc(k,1,nseq)
        y(1)=calfc(k,2,nseq)
        x(2)=calfc(m,1,nseq)
        y(2)=calfc(m,2,nseq)
        dist=SQRT((x(2)-x(1))**2+(y(2)-y(1))**2)
        if(dist.gt.ymax)then
        ymax=dist
            if(y(1).lt.y(2))then
            yx1=x(1)
            yx2=x(2)
            yy1=y(1)
            yy2=y(2)
            else
          yx1=x(2)
            yx2=x(1)
            yy1=y(2)
            yy2=y(1)
            endif
        endif
67          continue
68      continue

        write(16,*)'yx2: ',yx2,' yy2: ',yy2
            yx2=yx2-yx1
            yy2=yy2-yy1
        write(16,*)'yx2: ',yx2,' yy2: ',yy2
            
      d=SQRT((yx2)**2+(yy2)**2)
      sinfi=yx2/d
      cosfi=yy2/d
      write(16,*)'cosfi:',cosfi,' sinfi: ',sinfi
           
        yx2po=yx2*cosfi+yy2*(-sinfi)
        yy2po=yx2*sinfi+yy2*cosfi
      if(ABS(yx2po).lt.9.0E-6)yx2po=0.0
      if(ABS(yy2po).lt.9.0E-6)yy2po=0.0
      write(16,*)'turn around z axis'
      write(16,*)'yx2: ',yx2po,' yy2: ',yy2po

      do 70 k=natch(l-1)+1,natch(ll)
      xold=coord(k,1,nseq)-yx1
      yold=coord(k,2,nseq)-yy1
      coord(k,1,nseq)=xold*cosfi+yold*(-sinfi)
      coord(k,2,nseq)=xold*sinfi+yold*cosfi
      if(ABS(coord(k,1,nseq)).lt.9.0E-6)coord(k,1,nseq)=0.0
      if(ABS(coord(k,2,nseq)).lt.9.0E-6)coord(k,2,nseq)=0.0
70    continue

      do 73 k=ncach(l-1)+1,ncach(ll)
	do 72 n=1,bocz(k)
      xold=atomy(k,n,2)-yx1
      yold=atomy(k,n,3)-yy1
      atomy(k,n,2)=xold*cosfi+yold*(-sinfi)
      atomy(k,n,3)=xold*sinfi+yold*cosfi
      if(ABS(atomy(k,n,2)).lt.9.0E-6)atomy(k,n,2)=0.0
      if(ABS(atomy(k,n,3)).lt.9.0E-6)atomy(k,n,3)=0.0
72	continue
73    continue

      do 75 k=ncach(l-1)+1,ncach(ll)
      xold=calfc(k,1,nseq)-yx1
      yold=calfc(k,2,nseq)-yy1
      calfc(k,1,nseq)=xold*cosfi+yold*(-sinfi)
      calfc(k,2,nseq)=xold*sinfi+yold*cosfi
      if(ABS(calfc(k,1,nseq)).lt.9.0E-6)calfc(k,1,nseq)=0.0
      if(ABS(calfc(k,2,nseq)).lt.9.0E-6)calfc(k,2,nseq)=0.0
75    continue

	if(nonhet.eq.0)then
      do 77 k=1,nhet
      xold=hetco(k,1,nseq)-yx1
      yold=hetco(k,2,nseq)-yy1
      hetco(k,1,nseq)=xold*cosfi+yold*(-sinfi)
      hetco(k,2,nseq)=xold*sinfi+yold*cosfi
      if(ABS(hetco(k,1,nseq)).lt.9.0E-6)hetco(k,1,nseq)=0.0
      if(ABS(hetco(k,2,nseq)).lt.9.0E-6)hetco(k,2,nseq)=0.0
77    continue
	elseif(nonhet.ne.1)then
	write(16,850)'nonhet'
	endif
 
80      continue
      
      
      goto 100
          
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')

100   end

