SUBROUTINE readstr(nat,nres,nch,natch,ncach,coord,calfc,hetco,nhet,nseq,protein,atomy,bocz,linia1,linia2,filetp,elipsa,bezwodor,&
resnam,fresnam,amass,protch,nmr_ok,drop,nonhet,het_present,liganddef,ligand_id,ligand,nlig,readhet,nohet_ca,resid,plres_id,&
planecrd,planetype,planedef,nplane,numf,flin1,flin2,fatom,fbocz,fnres,fnat,fresid,chainid)
IMPLICIT none
!***********************   reads atoms coordinates ***************************

REAL(8)::calfc(15000,3,2),coord(100000,3,2),hetco(10000,3,2),atomy(15000,30,4),amass(15000,30),ligand(100,3),planecrd(4,3,2),&
fatom(15000,30,4),x,y,z
INTEGER::natch(0:50),ncach(0:50),nres,bocz(15000),chaintype,resid(15000),fresid(15000),nplane(4,2),fbocz(15000),fnres,em,filetp,&
nat,nch,nhet,nseq,elipsa,nmr_ok,liganddef,nlig,readhet,nohet_ca,planetype,planedef,numf,i,k,n,napl,nr,fnat,nonhet
character line*80,atomt*4,restype*3,linia1(15000,30)*30,linia2(15000,30)*24,resnam(15000)*3,fresnam(15000)*3,chid*3,ligand_id*3,&
plres_id*3,flin1(15000,30)*30,flin2(15000,30)*24,chainid(15000)*3
logical not_nres1,terok,bezwodor,protein,protch(50),drop,het_present,match,matchpl

napl=0
match=.FALSE.
matchpl=.FALSE.
not_nres1=.FALSE.
i=0
if(filetp.eq.1)then
  terok=.FALSE.
5	read(1,'(a)',end=100,err=99)line
!.wczytaj wspolrzedne atomow
	if(line(1:4).eq.'ATOM')then
		terok=.FALSE. 
        backspace 1
		read(1,800,err=98)atomt,restype,chid,nr,x,y,z
		if(bezwodor)then
			if(line(13:13).eq.'H'.OR.line(14:14).eq.'H'.OR.line(15:15).eq.'H'.OR.line(16:16).eq.'H')goto 5
		endif		
!	nie wczytuje podwojnych wpisow aminokwasow oznaczonych roznymi literami na polu 17
		if(line(17:17).eq.'B'.OR.line(17:17).eq.'2')goto 5
		nat=nat+1 
        coord(nat,1,nseq)=x
        coord(nat,2,nseq)=y
        coord(nat,3,nseq)=z

		if(atomt.eq.' N  ')then
      		chaintype=1
			protein=.TRUE.
			nres=nres+1
			resid(nres)=nr
			chainid(nres)=chid
			resnam(nres)=restype
		 	if(.NOT.not_nres1)goto 10
			bocz(nres-1)=i
10			not_nres1=.TRUE.
			i=4
			amass(nres,1)=14.0067
			atomy(nres,1,1)=2
        	atomy(nres,1,2)=x
        	atomy(nres,1,3)=y
        	atomy(nres,1,4)=z
			linia1(nres,1)=line(1:30)
			linia2(nres,1)=line(55:80)
		elseif(atomt.eq.' P  ')then
      		chaintype=2
        	nres=nres+1
			resnam(nres)=restype
			resid(nres)=nr
			chainid(nres)=chid
		 	if(.NOT.not_nres1)goto 15
			bocz(nres-1)=i
15			not_nres1=.TRUE.
			i=1
			amass(nres,2)=30.97376
			atomy(nres,1,1)=1
        	atomy(nres,1,2)=x
			atomy(nres,1,3)=y
        	atomy(nres,1,4)=z
			linia1(nres,1)=line(1:30)
			linia2(nres,1)=line(55:80)
        elseif(atomt.eq.' CA ')then
			if(.NOT.not_nres1)then
				nres=nres+1
				not_nres1=.TRUE.
			endif
       		calfc(nres,1,nseq)= coord(nat,1,nseq)
       		calfc(nres,2,nseq)= coord(nat,2,nseq)
       		calfc(nres,3,nseq)= coord(nat,3,nseq)
			amass(nres,2)=12.0107
			atomy(nres,2,1)=1
            atomy(nres,2,2)=x
            atomy(nres,2,3)=y
            atomy(nres,2,4)=z
			linia1(nres,2)=line(1:30)
	        linia2(nres,2)=line(55:80)
		elseif(atomt.eq.' C  ')then
			if(.NOT.not_nres1)then
				nres=nres+1
				not_nres1=.TRUE.
			endif
			amass(nres,3)=12.0107
			atomy(nres,3,1)=1
            atomy(nres,3,2)=x
            atomy(nres,3,3)=y
            atomy(nres,3,4)=z
			linia1(nres,3)=line(1:30)
            linia2(nres,3)=line(55:80)
		elseif(atomt.eq.' O  ')then
			if(.NOT.not_nres1)then
				nres=nres+1
				not_nres1=.TRUE.
			endif
			amass(nres,4)=15.9994
			atomy(nres,4,1)=3
		    atomy(nres,4,2)=x
            atomy(nres,4,3)=y
            atomy(nres,4,4)=z
			linia1(nres,4)=line(1:30)
            linia2(nres,4)=line(55:80)
		elseif(atomt.eq.' OXT')then
			if(.NOT.not_nres1)then
				nres=nres+1
				not_nres1=.TRUE.
			endif
			i=i+1
			amass(nres,4)=15.9994
			atomy(nres,i,1)=3
            atomy(nres,i,2)=x
            atomy(nres,i,3)=y
            atomy(nres,i,4)=z
			linia1(nres,i)=line(1:30)
            linia2(nres,i)=line(55:80)
        else
			if(.NOT.not_nres1)then
				nres=nres+1
				not_nres1=.TRUE.
			endif
            i=i+1
			atomy(nres,i,2)=x
            atomy(nres,i,3)=y
            atomy(nres,i,4)=z
			linia1(nres,i)=line(1:30)
            linia2(nres,i)=line(55:80)
!....check atom type, assigns atomic masses
			if(atomt(2:2).eq.'C')then
				amass(nres,i)=12.0107
				atomy(nres,i,1)=1
			elseif(atomt(2:2).eq.'N')then
				amass(nres,i)=14.0067
				atomy(nres,i,1)=2
			elseif(atomt(2:2).eq.'O')then
				amass(nres,i)=15.9994
				atomy(nres,i,1)=3
			elseif(atomt(2:2).eq.'S')then
				amass(nres,i)=32.065
				atomy(nres,i,1)=4
			elseif(atomt(2:2).eq.'H')then
				amass(nres,i)=1.00794
				atomy(nres,i,1)=5
			elseif(atomt(2:2).eq.'Ca')then
				amass(nres,i)=40.078
				atomy(nres,i,1)=6
			elseif(atomt(2:2).eq.'Na')then
				amass(nres,i)=22.989769
				atomy(nres,i,1)=7	
			elseif(atomt(2:2).eq.'K')then
				amass(nres,i)=39.0983
				atomy(nres,i,1)=8
			elseif(atomt(2:2).eq.'Cl')then
				amass(nres,i)=35.453
				atomy(nres,i,1)=9
			elseif(atomt(2:2).eq.'F')then
				amass(nres,i)=18.9984
				atomy(nres,i,1)=10
			elseif(atomt(2:2).eq.'P')then
				amass(nres,i)=30.97376
				atomy(nres,i,1)=11
			endif 
        endif
	 	goto 5
     endif
!...hetero atoms section
      
     if(line(1:6).eq.'HETATM')then
     	if(readhet.eq.1)then
        	backspace 1
     		read(1,800,err=98)atomt,restype,chid,nr,x,y,z
      		if(restype.eq.'HOH')goto 5
			het_present=.TRUE.
       		nhet=nhet+1
			if(liganddef.eq.1)then
			 	if(restype.eq.ligand_id)then
		 			nlig=nlig+1	
	 				match=.TRUE.
	 				ligand(nlig,1)=x
				 	ligand(nlig,2)=y
	 				ligand(nlig,3)=z
		 		endif
			elseif(liganddef.eq.0)then
		 		nlig=nlig+1	
		 		match=.TRUE.
		 		ligand(nlig,1)=x
		 		ligand(nlig,2)=y
		 		ligand(nlig,3)=z
			elseif(liganddef.ne.0)then
				write(16,850)'liganddef'
				stop
			endif!end of liganddef
	
			if(planedef.eq.1)then
				if(planetype.eq.1)then	 
					if(restype.eq.plres_id)then
						if(napl.lt.4)then
	 						napl=napl+1	
					 		matchpl=.TRUE.
					 		planecrd(napl,1,nseq)=x
							planecrd(napl,2,nseq)=y
					 		planecrd(napl,3,nseq)=z
				 		endif
	 				endif
				elseif(planetype.ne.0)then
					write(16,850)'planetype'
					stop
				endif
			elseif(planedef.ne.0)then
				write(16,850)'planedef'
				stop
			endif !end of planedef
			hetco(nhet,1,nseq)=x
			hetco(nhet,2,nseq)=y
			hetco(nhet,3,nseq)=z
			if(atomt.eq.' CA ')then
	 			if(nohet_ca.eq.0)then
        	    	nres=nres+1
       				calfc(nres,1,nseq)=hetco(nhet,1,nseq)
       				calfc(nres,2,nseq)=hetco(nhet,2,nseq)
       				calfc(nres,3,nseq)=hetco(nhet,3,nseq)
	 			elseif(nohet_ca.ne.1)then
	 				write(16,850)'nohet_ca'
	 			endif
    		goto 5
			endif!end of atomt
		endif!end of readhet
    endif!end of HETATM
    goto 5       

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...reads coordinates from x,y,z file    
elseif(filetp.eq.2)then
	k=1
80	read(1,803,end=100)x,y,z
	coord(k,1,nseq)=x
	coord(k,2,nseq)=y
	coord(k,3,nseq)=z
	k=k+1
	goto 80
elseif(filetp.eq.3)then
else
	write(16,850)'filetp'
	stop
endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!...end of reading section

98	write(16,*)'reading error',atomt,restype,chid,nr,x,y,z
	drop=.TRUE.
	return      
99    write(16,*)'reading err'
      drop=.TRUE.
      return
100   write(16,*)'coordinates read...OK'
	nch=nch-1
	
!...control
if(liganddef.eq.1)then
	if(.NOT.match)then
	write(16,*)'ligand in question: ',ligand_id, ' not found' 
	write(*,*)'ligand in question: ',ligand_id, ' not found' 
	endif
endif
	
if(planetype.eq.1)then
	if(.NOT.matchpl)then
	write(16,*)'residue in question: ',plres_id, ' not found' 
	write(*,*)'residue in question: ',plres_id, ' not found' 
	endif
elseif(planetype.ne.0)then
	write(*,*)planetype,'zupa2'
	write(16,850)'planetype'
	stop
endif 

!...defines atom #s based reference plane
if(planedef.eq.1)then
	if(planetype.eq.0)then
		do i=1,3
		planecrd(1,i,nseq)=coord(nplane(1,nseq),i,nseq)
		planecrd(2,i,nseq)=coord(nplane(2,nseq),i,nseq)
		planecrd(3,i,nseq)=coord(nplane(3,nseq),i,nseq)
		planecrd(4,i,nseq)=coord(nplane(4,nseq),i,nseq)			
		end do
	elseif(planetype.ne.1)then
		write(16,850)'planetype'
		stop
	endif
elseif(planedef.ne.0)then
	write(16,850)'planedef'
	stop
endif

!...finishing section, controls proper residue count in chains
		
nch=nch+1
bocz(nres)=i
i=0
ncach(nch)=nres
natch(nch)=nat
			
if(numf.eq.1)then
	fnres=nres
	fnat=nat
	do k=1,fnres
	fresid(k)=resid(k)
	fbocz(k)=bocz(k)
	fresnam(k)=resnam(k)
		do n=1,fbocz(k)
		fatom(k,n,1)=atomy(k,n,1)
		fatom(k,n,2)=atomy(k,n,2)
		fatom(k,n,3)=atomy(k,n,3)
		fatom(k,n,4)=atomy(k,n,4)
		flin1(k,n)=linia1(k,n)
		flin2(k,n)=linia2(k,n)
		end do	
	end do	
endif

write(16,*)'total number of residues: ',nres
write(16,*)'total number of atoms: ',nat
goto 150

800    format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
803   format(3f8.3)
850     format('Wrong value of following parameter: ',a8,/&
       'restart with proper value set')

150       end

