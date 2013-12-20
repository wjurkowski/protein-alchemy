SUBROUTINE vr(ws,nres)
IMPLICIT none
integer iq,nres,isp,i1,i2,i3,i4,i5,ip,il,ik,is,k,l,nat
real(8)::ws(15000,3),wr(15000,3),wa(15000,3),xc,yc,zc,xo,yo,zo,a,w,x11,x12,x13,x21,x22,x23,x31,x32,x33,&
x1,x2,x3,z1,z2,z3,rs,ra,dss,daa,ss,pr,xxx1,xxx2,x,y,z,wz,xa,xb,ya,yb,hb,ha,xx,yy,zz,dd,ddd,vm1,vm2,vp1,&
vp2,v01,v0,v11,v22,v0p1,v0m1,v0m1abs,v0p1abs,v0p2,v0m2,v0m2abs,v0p2abs,v11abs,v22abs,prl,vsr,vmp
!write(*,*)ws(1,1),ws(1,2),ws(1,3)
!write(*,*)ws(2,1),ws(2,2),ws(2,3)
!write(*,*)ws(3,1),ws(3,2),ws(3,3)
!write(*,*)ws(4,1),ws(4,2),ws(4,3)
!write(*,*)ws(5,1),ws(5,2),ws(5,3)

		nat=nres*5

     do 9 is=1,1
              isp=3
!	write(*,*)ws(32,1),ws(32,2),ws(32,3)
!	write(*,*)ws(37,1),ws(37,2),ws(37,3)
      do 13 iq=4,nat-20,5
!write(*,*)ws(isp,1)
              i1=iq+5
              i2=iq+10
              i3=iq+15
              i4=iq+20
              i5=iq
              xc=ws(i1,1)+ws(i2,1)+ws(i3,1)+ws(i4,1)+ws(i5,1)
              yc=ws(i1,2)+ws(i2,2)+ws(i3,2)+ws(i4,2)+ws(i5,2)
              zc=ws(i1,3)+ws(i2,3)+ws(i3,3)+ws(i4,3)+ws(i5,3)
              xc=xc/5.
              yc=yc/5.
              zc=zc/5.
              do 14 k=1,nat
              ws(k,1)=ws(k,1)-xc
              ws(k,2)=ws(k,2)-yc
              ws(k,3)=ws(k,3)-zc
14            continue
!	write(*,*)xc,yc,zc
!	write(*,*)ws(32,1),ws(32,2),ws(32,3)
!	write(*,*)ws(37,1),ws(37,2),ws(37,3)
!write(*,*)'dasdasdas'
              i1=iq+5+1
              i2=iq+10+1
              i3=iq+15+1
              i4=iq+20+1
              i5=iq+1
              xo=ws(i1,1)+ws(i2,1)+ws(i3,1)+ws(i4,1)+ws(i5,1)
              yo=ws(i1,2)+ws(i2,2)+ws(i3,2)+ws(i4,2)+ws(i5,2)
              zo=ws(i1,3)+ws(i2,3)+ws(i3,3)+ws(i4,3)+ws(i5,3)
              xo=xo/5.
              yo=yo/5.
              zo=zo/5.
              w=xo*xo+yo*yo
              a=w+zo*zo
              if(w.lt.0.00025) go to 3949
              a=sqrt(a)
              w=sqrt(w)
              x11=yo/w
              x12=zo*xo/(a*w)
              x13=xo/a
              x21=-xo/w
              x22=yo*zo/(w*a)
              x23=yo/a
              x31=0.0
              x32=-w/a
              x33=zo/a
              do 144 l=1,nat
              z1=ws(l,1)
              z2=ws(l,2)
              z3=ws(l,3)
              ws(l,1)=z1*x11+z2*x21+z3*x31
              ws(l,2)=z1*x12+z2*x22+z3*x32
              ws(l,3)=z1*x13+z2*x23+z3*x33
!	write(*,*)z1,z2,z3
144           continue
!	write(*,*)ws(32,1),ws(32,2),ws(32,3)
!	write(*,*)ws(37,1),ws(37,2),ws(37,3)
3949          continue
!write(*,*)iq,rs,dss,ss
              call sred (iq,rs,dss,ss,wr,nat,ws)
!write(*,*)iq,rs,dss,ss
              call anal (iq,ra,daa,a,wa,nat,ws)
!write(*,*)iq,ra,dss,daa
              if (dss.lt.daa) go to 400
              if (daa.lt.dss) go to 401
400           do 454 ip=1,nat
              do 434 il=1,3   
              ws(ip,il)=wr(ip,il)
434           continue 
454           continue
              pr=rs
              go to 466
401            continue
402            do 404 k=1,nat
               ws(k,1)=wa(k,1)
               ws(k,2)=wa(k,2)
               ws(k,3)=wa(k,3)
404            continue
               pr=ra
466             xxx1=0.0
                xxx2=0.0
               do 405 k=iq+9,iq+14               
                xxx1=xxx1+ws(k,1)/6.
                xxx2=xxx2+ws(k,2)/6.
405            continue
406            continue
               x=xxx1
               y=xxx2
               wz=x*x+y*y
!	write(*,*)wz
               if(wz.lt.0.00025) go to 4488
               wz=sqrt(wz)
               x11=x/wz
               x12=-y/wz
               x21=-x12
               x22=x11
               x33=1.0                    
               x31=0.0
               x32=0.0
               x13=0.0
               x23=0.0
               do 1471 l=1,nat
               x1=ws(l,1)
               x2=ws(l,2)
               x3=ws(l,3)
               ws(l,1)=x1*x11+x2*x21+x3*x31   
               ws(l,2)=x1*x12+x2*x22+x3*x32   
               ws(l,3)=x1*x13+x2*x23+x3*x33 
1471           continue
4488           continue
               xa=0.0
               xb=0.0
               ya=0.0
               yb=0.0         
               do 1492 ik=iq+4,iq+9
               xb=xb+ws(ik,1)/6. 
               yb=yb+ws(ik,2)/6.
1492           continue
               do 1493 ik=iq+14, iq+19
               xa=xa+ws(ik,1)/6.
               ya=ya+ws(ik,2)/6.
1493           continue
               hb=(yb*yb+xb*xb)
               hb=sqrt(hb)
               hb=1./hb
               hb=hb*xb
               ha=(ya*ya+xa*xa)
               ha=sqrt(ha)
               ha=1./ha
               ha=ha*xa
               hb=57.3*asin(hb)
               ha=57.3*asin(ha)
              x=ws(iq+10,1)
               y=ws(iq+10,2)
               z=ws(iq+10,3)
               xx=ws(iq+11,1)-x
               yy=ws(iq+11,2)-y
               zz=ws(iq+11,3)-z
               dd=xx*xx+yy*yy  
               ddd=dd+zz*zz
               dd=sqrt(dd)
               ddd=sqrt(ddd)
               v0=zz/ddd 
               v0=57.3*acos(v0) 
               x=ws(iq+5,1)
               y=ws(iq+5,2)
               z=ws(iq+5,3)
               xx=ws(iq+6,1)-x
               yy=ws(iq+6,2)-y
               zz=ws(iq+6,3)-z
               dd=xx*xx+yy*yy  
               ddd=dd+zz*zz
               dd=sqrt(dd)
               ddd=sqrt(ddd)
               vm1=zz/ddd 
               vm1=57.3*acos(vm1)
               x=ws(iq,1)
               y=ws(iq,2)
               z=ws(iq,3)
               xx=ws(iq+1,1)-x
               yy=ws(iq+1,2)-y
               zz=ws(iq+1,3)-z
               dd=xx*xx+yy*yy  
               ddd=dd+zz*zz
               dd=sqrt(dd)
               ddd=sqrt(ddd)
               vm2=zz/ddd 
               vm2=57.3*acos(vm2) 
               x=ws(iq+15,1)
               y=ws(iq+15,2)
               z=ws(iq+15,3)
               xx=ws(iq+16,1)-x
               yy=ws(iq+16,2)-y
               zz=ws(iq+16,3)-z
               dd=xx*xx+yy*yy  
               ddd=dd+zz*zz
               dd=sqrt(dd)
               ddd=sqrt(ddd)
               vp1=zz/ddd 
               vp1=57.3*acos(vp1)
               x=ws(iq+20,1)
               y=ws(iq+20,2)
               z=ws(iq+20,3)
               xx=ws(iq+21,1)-x
               yy=ws(iq+21,2)-y
               zz=ws(iq+21,3)-z
               dd=xx*xx+yy*yy  
               ddd=dd+zz*zz
               dd=sqrt(dd)
               ddd=sqrt(ddd)
               vp2=zz/ddd 
               vp2=57.3*acos(vp2)
300            format(5f10.3)
               v22=vm2-vp2
               v22abs=abs(vm2-vp2)
               v11=vm1-vp1
               v11abs=abs(vm1-vp1)
               v0m1=v0-vm1
               v0m1abs=abs(v0-vm1)
               v0p1=v0-vp1
               v0p1abs=abs(v0-vp1)
               v0m2=v0-vm2
               v0m2abs=abs(v0-vm2)
               v0p2=v0-vp2
               v0p2abs=abs(v0-vp2)
               prl=log(pr)
               v01=abs(v0)-abs(vp1)   
               vmp=abs(vm1)-abs(vp1)
               vsr=(abs(vm1)+abs(vp1))/2.
!write(*,*)'is:',is
!write(*,*)isp,pr,prl,vsr
if(is.eq.1)then 
!	write(*,*)'dupablada'
	write(4,498) isp,pr,prl,vsr
endif
     
       if(vsr.gt.180.) write (6,*) vm1,vp1,v0
               isp=isp+1
13    continue
499            format(5f10.3)
498            format(i5,4f16.3)
               rewind 12
9    continue

RETURN

END                
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE sred (iq,rs,dss,ss,wr,nat,ws)
real(8)::p(2),ws(15000,3),wr(15000,3),ss,dss,xx,yy,rs,rr,x,y,wz,x11,x12,x13,x21,x22,x23,x31,x32,x33,&
z1,z2,z3
               do 7 i=1,nat
               do 6 ii=1,3
               wr(i,ii)=ws(i,ii)
6              continue
!	write(*,*)wr(i,1),wr(i,2),wr(i,3)
7              continue
!	write(*,*)wr(1,1),wr(1,2),wr(1,3)
               do 8 ii=1,2
               p(ii)=0.0
 8             continue
               ip=iq-1
               ik=iq+24
               do 11 il=1,2
               do 12 im=ip,ik
               p(il)=p(il)+wr(im,il)/26.
!	write(*,*)p(il)
12             continue
11             continue
!	write(*,*)p(1),p(2)
               do 14 i=1,nat
               do 15 im=1,2     
               wr(i,im)=wr(i,im)-p(im)
15             continue
14             continue
!	write(*,*)wr(1,1),wr(1,2)
               ss=0.0
               dss=0.0
               do 16 i=ip,ik
               xx=wr(i,1)
               yy=wr(i,2)
               ss=ss+xx*xx+yy*yy 
!	write(*,*)ss
16             continue
               rs=ss/26.
               rs=sqrt(rs)
               do 17 i=ip,ik
               xx=wr(i,1)
               yy=wr(i,2)
               rr=xx*xx+yy*yy
               rr=sqrt(rr)
               dss=dss+(rr-rs)*(rr-rs)
17             continue
               ip=iq+9
               ik=iq+14
                x=0.0
                y=0.0
               do 99 i=ip,ik
               x=x+wr(i,1)/6.
               y=y+wr(i,2)/6.
99             continue
               wz=x*x+y*y
               if(wz.lt.0.00025) go to 1372
                wz=sqrt(wz)
                x11=x/wz
                x12=-y/wz
                x13=0.0
                x21=-x12
                x22=x11
                x23=0.0
                x31=0.0
                x32=0.0
                x33=1.0
                do 1371 l=1,nat
                z1=wr(l,1)
                z2=wr(l,2)
                z3=wr(l,3)
                wr(l,1)=z1*x11+z2*x21+z3*x31
                wr(l,2)=z1*x12+z2*x22+z3*x32
                wr(l,3)=z1*x13+z2*x23+z3*x33
1371            continue
1372            return 
               end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE anal (iq,ra,daa,a,wa,nat,ws)
real(8)::p(2),wa(15000,3),pp(5,2),ws(15000,3),wr(15000,3),ss,daa,xx,yy,ra,rr,x,y,wz,x11,x12,x13,x21,&
x22,x23,x31,x32,x33,x1,x2,x3,y1,y2,y3,z1,z2,z3,g1,g2,bb,ssxx,aa,tt
integer nat

               do 6 i=1,nat
               do 5 ij=1,3
               wa(i,ij)=ws(i,ij)
5              continue
6              continue
               do 8 ii=1,2
               p(ii)=0.0
               do 44 ik=1,5
               pp(ik,ii)=0.0
44             continue
 8             continue
               ip=iq-3
               ik=iq+21
               do 11 il=1,2
               do 12 im=ip,ik
               p(il)=p(il)+wa(im,il)/25.
12             continue
11             continue
               do 14 i=1,nat
               do 15 im=1,2     
               wa(i,im)=wa(i,im)-p(im)
15             continue
14             continue
               ip=iq-1
               ik=iq+4
               do 100 iz=1,3
               do 120 im=ip,ik
               pp(iz,1)=pp(iz,1)+wa(im,1)/6.
               pp(iz,2)=pp(iz,2)+wa(im,2)/6.
120             continue
               ip=ip+10
               ik=ik+10
100             continue
               x1=pp(1,1)
               x2=pp(2,1)
               x3=pp(3,1)
               y1=pp(1,2)
               y2=pp(2,2)
               y3=pp(3,2)
               g1=(x1-x3)*(x1*x1-x2*x2)-(x1*x1-x3*x3)*(x1-x2)
               g2=(x1-x3)*(y1*y1-y2*y2)-(y1*y1-y3*y3)*(x1-x2)
               bb=2.*((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))
               bb=1./bb
               bb=(g1+g2)*bb
               ssxx=(x1-x2)
               ssxx=abs(ssxx)
               aa=2.*(x1-x2)
               aa=1./aa
               aa=aa*((x1*x1-x2*x2)+(y1*y1-y2*y2)-2.*bb*(y1-y2))
               do 414 id=1,nat
               wa(id,1)=wa(id,1)-aa
               wa(id,2)=wa(id,2)-bb
414            continue
               do 771 i0=1,5
               do 772 j0=1,2
               pp(i0,j0)=0.0
772            continue
771            continue
               ip=iq-1
               ik=iq+4
               do 144 iz=1,3
               do 142 im=ip,ik
               pp(iz,1)=pp(iz,1)+wa(im,1)/6.
               pp(iz,2)=pp(iz,2)+wa(im,2)/6.
142            continue
               ip=ip+10
               ik=ik+10
144            continue
               ra=pp(1,1)*pp(1,1)+pp(1,2)*pp(1,2)
               ra=sqrt(ra)
               daa=0.0
               ip=iq-1
               ik=iq+21
               do 16 i=ip,ik
               xx=wa(i,1)
               yy=wa(i,2)
               tt=xx*xx+yy*yy
               tt=sqrt(tt) 
               daa=daa+(tt-ra)*(tt-ra) 
661            format(1x,'wartosci tt i ra',3f10.3)
16             continue
               ip=iq+9
               ik=iq+14
               x=0.0
               y=0.0
               do 99 i=ip,ik
               x=x+wa(i,1)/6.
               y=y+wa(i,2)/6.
99             continue
               wz=x*x+y*y
               if(wz.lt.0.00025) go to 1372
                wz=sqrt(wz)
                x11=x/wz
                x12=-y/wz
                x13=0.0
                x21=-x12
                x22=x11
                x23=0.0
                x31=0.0
                x32=0.0
                x33=1.0
                do 1371 l=1,nat
                z1=wa(l,1)
                z2=wa(l,2)
                z3=wa(l,3)
                wa(l,1)=z1*x11+z2*x21+z3*x31
                wa(l,2)=z1*x12+z2*x22+z3*x32
                wa(l,3)=z1*x13+z2*x23+z3*x33
1371            continue
1372            return 
END
