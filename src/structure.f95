SUBROUTINE structure(atomy,bocz,atpo,resid,ncach,nch,minres,ellcrd,nres_str,edihass,estrtype,aamax,resnam_str,&
resnam_aa,resid_aa,fipsi,ellfipsi,namcor)
USE tablice
IMPLICIT none
	
REAL(8)::ellfipsi(15000,2),atomy(15000,30,4),fipsi(15000,3),atpo(15000,30,4),cooN(3),cooCA(3),cooC(3),&
 cooO(3),cooN2(3),aamax(20,14),cooCA2(3),cooC2(3),cooO2(3),cooN3(3),wxy_ck,wxz_ck,sina,cosa,x,y,z,&
shift1(3),shift2(3),ekat,ekatrad,sinell,cosell,sinb,cosb,prefi,prepsi,wxy_ca,wxz_ca,transx,transy,&
transz
INTEGER::bocz(15000),ncach(0:50),resid(15000),nres_str(20),resid_aa(20,1000),nch,minres,ellcrd,edihass,estrtype,&
l,m,n,ii,kpos,l1,l2
character resnam_str(20,1000)*1,resnam_aa(20,1000)*3,namcor*50	

!...creates structures either from scratch or from the native structures	
	
		
!0000000000000000000000000000000000000000000000000000000000000000000000
!...creates structures from scratch	
if(estrtype.eq.2)then
write(16,*)'writes elliptic structure based on: read in structural sequence...'

	do 65 i=1,1
!=====================================================================
!...fi,psi dihedrals assignment basing on structural sequence
	if(edihass.eq.1)then
!...ads eight residues more for motifes flanking
	nres_str(i)=nres_str(i)+8

	do 1 k=1,4
	ellfipsi(k,1)=180.000
	ellfipsi(k,2)=180.000
1	continue
	do 4 k=5,nres_str(i)-4
	 do 3 m=1,20
	 if(resnam_aa(i,k-4).eq.aminok(m))then	
	  do 2 ii=1,7
	  if(resnam_str(i,k-4).eq.AS(ii))then
	  l1=(ii-1)*2+1
	  l2=(ii-1)*2+2	 
	  ellfipsi(k,1)=aamax(m,l1)
	  ellfipsi(k,2)=aamax(m,l2)
	write(*,*)'aaa',ellfipsi(k,1),ellfipsi(k,2)
	  goto 4	
	  endif
2	  continue
	 endif	
3	 continue	
4	continue
	do 6 k=nres_str(i)-3,nres_str(i)
	ellfipsi(k,1)=180.000
	ellfipsi(k,2)=180.000
6	continue
	elseif(edihass.ne.0)then
      write(16,850)'edihass'
      stop
      endif
	
!====================================================================	
!...creates structure from scratch - predefines backbone coordinates
	
	cooN(1)=0.039
	cooN(2)=-0.028
	cooN(3)=0.000
	cooCA(1)=1.499
	cooCA(2)=-0.043
	cooCA(3)=0.000
	cooC(1)=2.055
	cooC(2)=1.361
	cooC(3)=0.000
	cooO(1)=1.321		
	cooO(2)=2.356
	cooO(3)=0.011
	cooN2(1)=3.371			 
	cooN2(2)=1.462
	cooN2(3)=0.000
	cooCA2(1)=4.047		
	cooCA2(2)=2.756
	cooCA2(3)=0.000
	cooC2(1)=5.547		
	cooC2(2)=2.582
	cooC2(3)=0.000
	cooO2(1)=6.078		
	cooO2(2)=1.465
	cooO2(3)=-0.011
	cooN3(1)=6.259			 
	cooN3(2)=3.694
	cooN3(3)=0.000
		
	cooN(1)=0.039
	cooN(2)=-0.028
	cooN(3)=0.000
	cooCA(1)=1.499
	cooCA(2)=-0.043
	cooCA(3)=0.000
	cooC(1)=2.055
	cooC(2)=1.361
	cooC(3)=0.000
	cooO(1)=1.321		
	cooO(2)=2.356
	cooO(3)=0.011
	cooN2(1)=3.371			 
	cooN2(2)=1.462
	cooN2(3)=0.000
	cooCA2(1)=4.047		
	cooCA2(2)=2.756
	cooCA2(3)=0.000
	cooC2(1)=5.547		
	cooC2(2)=2.582
	cooC2(3)=0.000
	cooO2(1)=6.078		
	cooO2(2)=1.465
	cooO2(3)=-0.011
	
	atpo(1,1,2)=cooN(1)
	atpo(1,1,3)=cooN(2)
	atpo(1,1,4)=cooN(3)
	atpo(1,2,2)=cooCA(1)
	atpo(1,2,3)=cooCA(2)
	atpo(1,2,4)=cooCA(3)
	atpo(1,3,2)=cooC(1)
	atpo(1,3,3)=cooC(2)
	atpo(1,3,4)=cooC(3)
	atpo(1,4,2)=cooO(1)
	atpo(1,4,3)=cooO(2)
	atpo(1,4,4)=cooO(3)
	atpo(2,1,2)=cooN2(1)
	atpo(2,1,3)=cooN2(2)
	atpo(2,1,4)=cooN2(3)
	atpo(2,2,2)=cooCA2(1)
	atpo(2,2,3)=cooCA2(2)
	atpo(2,2,4)=cooCA2(3)
	atpo(2,3,2)=cooC2(1)
	atpo(2,3,3)=cooC2(2)
	atpo(2,3,4)=cooC2(3)
	atpo(2,4,2)=cooO2(1)
	atpo(2,4,3)=cooO2(2)
	atpo(2,4,4)=cooO2(3)
	transx=cooN3(1)-cooN(1)
	transy=cooN3(2)-cooN(2)
	transz=cooN3(3)-cooN(3)
	prefi=180.000
	prepsi=180.000

	do 10 k=3,nres_str(i),2
	 do 9 n=1,4
	atpo(k,n,2)=atpo(1,n,2)+(((k-1)/2)*transx)
	atpo(k,n,3)=atpo(1,n,3)+(((k-1)/2)*transy)
	atpo(k,n,4)=atpo(1,n,4)+(((k-1)/2)*transz)
	atpo(k+1,n,2)=atpo(2,n,2)+(((k-1)/2)*transx)
	atpo(k+1,n,3)=atpo(2,n,3)+(((k-1)/2)*transy)
	atpo(k+1,n,4)=atpo(2,n,4)+(((k-1)/2)*transz)
9	 continue	
10	continue
	
!====================================================================
!...eliptic structure creation

!...transformes consecutive residues
	do 60 k=1,nres_str(i)-1
!...psi
	shift1(1)=atpo(k,2,2)
	shift1(2)=atpo(k,2,3)
	shift1(3)=atpo(k,2,4)
	

!.....shift and turn
	do 32 l=1,nres_str(i)
	do 31 n=1,4
	atpo(l,n,2)=atpo(l,n,2)-shift1(1)
	atpo(l,n,3)=atpo(l,n,3)-shift1(2)
	atpo(l,n,4)=atpo(l,n,4)-shift1(3)
31	continue
32	continue

!.......rzut na os x	
	wxy_ck=SQRT(atpo(k,3,2)**2+atpo(k,3,3)**2)
	sina=atpo(k,3,3)/wxy_ck
	cosa=atpo(k,3,2)/wxy_ck
	if(ABS(sina).lt.9.0E-7)sina=0.0
	if(ABS(cosa).lt.9.0E-7)cosa=0.0

	do 34 l=1,nres_str(i)
	do 33 n=1,4
	x=atpo(l,n,2)*cosa+atpo(l,n,3)*sina
	y=-atpo(l,n,2)*sina+atpo(l,n,3)*cosa
	if(ABS(y).lt.9.0E-7)y=0.0
	atpo(l,n,2)=x
	atpo(l,n,3)=y
33      continue
34      continue

!.......rzut na os z
	wxz_ck=SQRT(atpo(k,3,2)**2+atpo(k,3,4)**2)
	sinb=atpo(k,3,2)/wxz_ck
	cosb=atpo(k,3,4)/wxz_ck

	do 36 l=1,nres_str(i)
        do 35 n=1,4	
	x=atpo(l,n,2)*cosb-atpo(l,n,4)*sinb
	z=atpo(l,n,2)*sinb+atpo(l,n,4)*cosb
	if(ABS(x).lt.9.0E-7)x=0.0
	atpo(l,n,2)=x
	atpo(l,n,4)=z
35      continue
36      continue

!......ustalenie kata obrotu
	ekat=ellfipsi(k,2)-prepsi
	ekatrad=ekat*(PI/180)
	sinell=-SIN(ekatrad)
	cosell=-COS(ekatrad)
	
!......obrot do struktury eliptycznej
!.........wegiel i tlen karbonylowy k-tej reszty
	
	x=atpo(k,3,2)*cosell+atpo(k,3,3)*sinell
	y=-atpo(k,3,2)*sinell+atpo(k,3,3)*cosell
	atpo(k,3,2)=x
	atpo(k,3,3)=y
	x=atpo(k,4,2)*cosell+atpo(k,4,3)*sinell
	y=-atpo(k,4,2)*sinell+atpo(k,4,3)*cosell
	atpo(k,4,2)=x
	atpo(k,4,3)=y

!.........wszystkie atomy reszty k+1  	
	do 45 l=(k+1),nres_str(i)  
	do 44 n=1,4
	x=atpo(l,n,2)*cosell+atpo(l,n,3)*sinell
	y=-atpo(l,n,2)*sinell+atpo(l,n,3)*cosell
	atpo(l,n,2)=x
	atpo(l,n,3)=y
44	continue
45	continue

!...fi

	shift2(1)=atpo((k+1),1,2)
	shift2(2)=atpo((k+1),1,3)
	shift2(3)=atpo((k+1),1,4)

!.....shift and turn
	do 52 l=1,nres_str(i)
	do 51 n=1,4
	atpo(l,n,2)=atpo(l,n,2)-shift2(1)
        atpo(l,n,3)=atpo(l,n,3)-shift2(2)
        atpo(l,n,4)=atpo(l,n,4)-shift2(3)
51	continue
52	continue

!.......rzut na os x
	wxy_ca=SQRT(atpo(k+1,2,2)**2+atpo(k+1,2,3)**2)
	sina=atpo(k+1,2,3)/wxy_ca
	cosa=atpo(k+1,2,2)/wxy_ca
	
	do 54 l=1,nres_str(i)
	do 53 n=1,4
	x=atpo(l,n,2)*cosa+atpo(l,n,3)*sina
        y=-atpo(l,n,2)*sina+atpo(l,n,3)*cosa
        if(ABS(y).lt.9.0E-7)y=0.0
        atpo(l,n,2)=x
        atpo(l,n,3)=y
53	continue
54	continue

!.......rzut na os z
	wxz_ca=SQRT(atpo(k+1,2,2)**2+atpo(k+1,2,4)**2)
	sinb=atpo(k+1,2,2)/wxz_ca
	cosb=atpo(k+1,2,4)/wxz_ca
	
	do 56 l=1,nres_str(i)
	do 55 n=1,4
	x=atpo(l,n,2)*cosb-atpo(l,n,4)*sinb
	z=atpo(l,n,2)*sinb+atpo(l,n,4)*cosb
	if(ABS(x).lt.9.0E-7)x=0.0
	atpo(l,n,2)=x
	atpo(l,n,4)=z
55	continue
56	continue
	
!......ustalenie kata obrotu
	ekat=ellfipsi(k+1,1)-prefi
	ekatrad=ekat*(PI/180)
	sinell=-SIN(ekatrad)
	cosell=-COS(ekatrad)
	
!......obrot do struktury eliptycznej
!.........od c-beta reszty k+1
	do 57 n=3,4 	
	x=atpo(k+1,n,2)*cosell+atpo(k+1,n,3)*sinell
	y=-atpo(k+1,n,2)*sinell+atpo(k+1,n,3)*cosell
	atpo(k+1,n,2)=x
	atpo(k+1,n,3)=y
57	continue

!.........wszystkie atomy od reszty k+2  	
	do 59 l=(k+2),nres_str(i)  
	 do 58 n=1,4
	x=atpo(l,n,2)*cosell+atpo(l,n,3)*sinell
	y=-atpo(l,n,2)*sinell+atpo(l,n,3)*cosell
	atpo(l,n,2)=x
	atpo(l,n,3)=y
58	 continue
59	continue
		
60	 continue	
65	continue
!============================================================

!0000000000000000000000000000000000000000000000000000000000000000000000
!...creates structures from native structures
elseif(estrtype.eq.1)then
!=======================================================================
!...eliptic structure creation basing on the native structure coordinates
	do 80 m=1,nch
	  if((ncach(m)-ncach(m-1)).lt.minres)goto 80

	do 79 k=ncach(m-1)+1,ncach(m)
	 do 78 n=1,bocz(k)
	atpo(k,n,2)=atomy(k,n,2)
	atpo(k,n,3)=atomy(k,n,3)
	atpo(k,n,4)=atomy(k,n,4)
78	 continue
79	continue
80	continue

!...transforming consecutive residues
	do 125 m=1,nch
	  if((ncach(m)-ncach(m-1)).lt.minres)goto 125

	do 120 k=ncach(m-1)+1,ncach(m)-1
!...psi
	shift1(1)=atpo(k,2,2)
	shift1(2)=atpo(k,2,3)
	shift1(3)=atpo(k,2,4)
	

!.....shift and turn
	do 92 l=ncach(m-1)+1,ncach(m)
	do 91 n=1,bocz(l)
	atpo(l,n,2)=atpo(l,n,2)-shift1(1)
	atpo(l,n,3)=atpo(l,n,3)-shift1(2)
	atpo(l,n,4)=atpo(l,n,4)-shift1(3)
91	continue
92	continue

!.......rzut na os x	
	wxy_ck=SQRT(atpo(k,3,2)**2+atpo(k,3,3)**2)
	sina=atpo(k,3,3)/wxy_ck
	cosa=atpo(k,3,2)/wxy_ck
	if(ABS(sina).lt.9.0E-7)sina=0.0
	if(ABS(cosa).lt.9.0E-7)cosa=0.0

	do 94 l=ncach(m-1)+1,ncach(m)
	do 93 n=1,bocz(l)
	x=atpo(l,n,2)*cosa+atpo(l,n,3)*sina
	y=-atpo(l,n,2)*sina+atpo(l,n,3)*cosa
	if(ABS(y).lt.9.0E-7)y=0.0
	atpo(l,n,2)=x
	atpo(l,n,3)=y
93      continue
94      continue

!.......rzut na os z
	wxz_ck=SQRT(atpo(k,3,2)**2+atpo(k,3,4)**2)
	sinb=atpo(k,3,2)/wxz_ck
	cosb=atpo(k,3,4)/wxz_ck

	do 96 l=ncach(m-1)+1,ncach(m)
        do 95 n=1,bocz(l)	
	x=atpo(l,n,2)*cosb-atpo(l,n,4)*sinb
	z=atpo(l,n,2)*sinb+atpo(l,n,4)*cosb
	if(ABS(x).lt.9.0E-7)x=0.0
	atpo(l,n,2)=x
	atpo(l,n,4)=z
95      continue
96      continue

!......ustalenie kata obrotu
	ekat=ellfipsi(k,2)-fipsi(k,2)
	ekatrad=ekat*(PI/180)
	sinell=-SIN(ekatrad)
	cosell=-COS(ekatrad)
	
!......obrot do struktury eliptycznej
!.........wegiel i tlen karbonylowy k-tej reszty
	
	x=atpo(k,3,2)*cosell+atpo(k,3,3)*sinell
	y=-atpo(k,3,2)*sinell+atpo(k,3,3)*cosell
	atpo(k,3,2)=x
	atpo(k,3,3)=y
	x=atpo(k,4,2)*cosell+atpo(k,4,3)*sinell
	y=-atpo(k,4,2)*sinell+atpo(k,4,3)*cosell
	atpo(k,4,2)=x
	atpo(k,4,3)=y

!.........wszystkie atomy reszty k+1  	
	do 105 l=(k+1),ncach(m)  
	do 104 n=1,bocz(l)
	x=atpo(l,n,2)*cosell+atpo(l,n,3)*sinell
	y=-atpo(l,n,2)*sinell+atpo(l,n,3)*cosell
	atpo(l,n,2)=x
	atpo(l,n,3)=y
104	continue
105	continue

!...fi
	shift2(1)=atpo((k+1),1,2)
	shift2(2)=atpo((k+1),1,3)
	shift2(3)=atpo((k+1),1,4)

!.....shift and turn
	do 112 l=ncach(m-1)+1,ncach(m)
	do 111 n=1,bocz(l)
	atpo(l,n,2)=atpo(l,n,2)-shift2(1)
      atpo(l,n,3)=atpo(l,n,3)-shift2(2)
      atpo(l,n,4)=atpo(l,n,4)-shift2(3)
111	continue
112	continue

!.......rzut na os x
	wxy_ca=SQRT(atpo(k+1,2,2)**2+atpo(k+1,2,3)**2)
	sina=atpo(k+1,2,3)/wxy_ca
	cosa=atpo(k+1,2,2)/wxy_ca
	
	do 114 l=ncach(m-1)+1,ncach(m)
	do 113 n=1,bocz(l)
	x=atpo(l,n,2)*cosa+atpo(l,n,3)*sina
        y=-atpo(l,n,2)*sina+atpo(l,n,3)*cosa
        if(ABS(y).lt.9.0E-7)y=0.0
        atpo(l,n,2)=x
        atpo(l,n,3)=y
113	continue
114	continue

!.......rzut na os z
	wxz_ca=SQRT(atpo(k+1,2,2)**2+atpo(k+1,2,4)**2)
	sinb=atpo(k+1,2,2)/wxz_ca
	cosb=atpo(k+1,2,4)/wxz_ca
	
	do 116 l=ncach(m-1)+1,ncach(m)
	do 115 n=1,bocz(l)
	x=atpo(l,n,2)*cosb-atpo(l,n,4)*sinb
	z=atpo(l,n,2)*sinb+atpo(l,n,4)*cosb
	if(ABS(x).lt.9.0E-7)x=0.0
	atpo(l,n,2)=x
	atpo(l,n,4)=z
115	continue
116	continue
	
!......ustalenie kata obrotu
	ekat=ellfipsi(k+1,1)-fipsi(k+1,1)
	ekatrad=ekat*(PI/180)
	sinell=-SIN(ekatrad)
	cosell=-COS(ekatrad)
	
!......obrot do struktury eliptycznej
!.........od c-beta reszty k+1
	do 117 n=3,bocz(k+1) 	
	x=atpo(k+1,n,2)*cosell+atpo(k+1,n,3)*sinell
	y=-atpo(k+1,n,2)*sinell+atpo(k+1,n,3)*cosell
	atpo(k+1,n,2)=x
	atpo(k+1,n,3)=y
117	continue

!.........wszystkie atomy od reszty k+2  	
	do 119 l=(k+2),ncach(m)  
	 do 118 n=1,bocz(l)
	x=atpo(l,n,2)*cosell+atpo(l,n,3)*sinell
	y=-atpo(l,n,2)*sinell+atpo(l,n,3)*cosell
	atpo(l,n,2)=x
	atpo(l,n,3)=y
118	 continue
119	continue
		
120	 continue	
125	continue

!...saves new elliptic coordinates
  if(ellcrd.eq.1)then
	do 132 m=1,nch
	  if((ncach(m)-ncach(m-1)).lt.minres)goto 132
      open(18,file=namcor(1:(kpos-1))//"."//sfx(m)//'.ellcrd',status='unknown')
	  do k=ncach(m-1)+1,ncach(m)
		do n=1,bocz(k)
		  write(18,801)atpo(k,n,2),atpo(k,n,3),atpo(k,n,4)
		end do
	  end do
	  close(18)
132	continue
	write(16,*)'elliptic from native structure creation...OK'
  elseif(ellcrd.ne.0)then
	write(16,850)'ellcrd'
	stop
  endif
else
  write(*,850)'estrtype'
  stop
endif

goto 200
	
190	write(16,*)'ERR opening unit 3'	
801	format(3f8.3)	
803   format(a4,1x,i6,2x,a2,2x,a3,i2,1x,i3,4x,3f8.3,3x,i1,20x)
850	format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')
     	
200	end
