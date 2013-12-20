PROGRAM pr_alchem
!************************  protein structure analyser  ***********************
!........ Wiktor Jurkowski     
USE tablice
IMPLICIT none
character*50 ans*1,infile,fnamcor,namcor,wyniki*15,paramet*15,output*15,linia1(15000,30)*30,&
linia2(15000,30)*24,resnam(15000)*3,fresnam(15000)*3,faaseq(15000)*1,aaseq(15000)*1,asseq(15000)*1,&
flin1(15000,30)*30,flin2(15000,30)*24,inpalign,resnam_str(20,1000)*1,resnam_aa(20,1000)*3,aasekw,&
assekw,ligand_id*3,plres_id*3,data*8,czas*10,strefa*5,BUFFER*100,pdb1*100,pdb2,chainid(15000)*3
     
logical fpresent,fpresent2,refine,bezwodor,protein,protch(50),drop,het_present
          
integer natch(0:50),ncach(0:50),bocz(15000),ngrcv(-180:180,-180:180,20),nmot,resid(15000),resmap(15000),&
fbocz(15000),agrid(-30:30,-30:30,-30:30,2),num_mol,natt,fresid(15000),gpk(-30:30,-30:30,-30:30,2),gpn(-30:30,-30:30,-30:30,2),&
fnres,resid_str(20,1000),resid_aa(20,1000),nres_aa(20),nres_str(20),nrestab(100),anal_bs,bsalign,cmap,cmaptype,dccaprof,fnat,&
dihsave,edihass,elipsa,elldih,ellcrd,ellpdb,eceppok,early_vol,estrtype,filetp,final,fold,gridok,gridtype,kcenter,dorama,rtype,&
kwsize,kwsize2,liganddef,mostki,motdl,motifpdb,minres,ngcat,wart(8),nmr_ok,nohydr,nonhet,nohet_ca,outpdb,overlay,pdbsave,planedef,&
planetype,readhet,residue_ref,refsys,savefit,dihstat,stralf,stranal,trnsys,torsje,vear,volfit,volref,nplane(4,2),kwezel(3),kwezelx,&
kwezely,kwezelz,j,l,m,n,kpos,kpos2,lpoz2,lpoz5,lp,lpca,numf,nat,algn_typ,AAoverlay,tekrok,nhet,nasseq,nch,nlig,nres,nr,nseq,meanr,&
nseqaln,molbm(1000),saveseq,seqtype,seq_motif2(100,20,2),bcutoff(20),ecutoff(20),alnall,aadihsave,ignmtf,block,bbaln,cmapmol,&
inpseqtype,fit,dih_outliers,capri,recs,rece,ligs,lige
     
real(8)::coord(100000,3,2),calfc(15000,3,2),hetco(10000,3,2),am(15000,30),atomy(15000,30,4),fipsi(15000,3),&
ellfipsi(15000,2),ligand(100,3),tstart,tstop,time,ellatomy(15000,30,4),bscenter(3),planecrd(4,3,2),fatom(15000,30,4),seed(3),&
grdtab(15000,30),grdcon(15000,30),wind2(100,10000),tfrepl(15000,30),cutoff_bs,cmapcut,emesh,grmesh,spacex,spacey,spacez,&
ucenterx,ucentery,ucenterz,space(3),ucenter(3),aamax(20,14),ws(15000,3),d,x1,x2,y1,y2,z1,z2

!00000000000000000000000000000000000000000000000000000000000000000000
!..starting section
call DATE_AND_TIME(data,czas,strefa,wart)
call CPU_TIME(tstart)

open(16,file='main.out',access='APPEND',status='unknown',action='readwrite',err=100)
write(16,800)
800     format(' '/&
        '	########  SIATAMI - Protein structure analyser ########'/&
        ' '/)
write(16,'(A25,1X,I4,A1,I2,A1,I2,1X,A2,1X,I2,A1,I2,1X,A10)')'Started interactively on:',wart(1),'-',&
wart(2),'-',wart(3),'at',wart(5),':',wart(6),'local time'
write(16,*)' '
write(16,*)'Starting options used:' 
     	

!00000000000000000000000000000000000000000000000000000000000000000000000
!..program flow parameters
!....controls for parameters
CALL GETARG(1,BUFFER)
write(16,*)'Parameter file used: ',BUFFER
inquire(file=BUFFER,exist=fpresent)
if(.not.fpresent)then
	write(*,*)'Program use:'
	write(*,*)'pr_alchem1 param pdb1 [pdb2]'
	stop
endif
READ(BUFFER,*)paramet

write(16,*)'File with program flow control parameters: ',paramet
write(16,*)'reads in parameters...'
!...reads in parameters for analysis of sequences and structures	
open(2,file=paramet,status='old',err=100) 
call readparam(filetp,nohydr,trnsys,refsys,nplane,torsje,elipsa,elldih,tekrok,volfit,volref,kcenter,gridok,ngcat,kwezel,space,&
ucenter,pdbsave,nonhet,mostki,minres,nmr_ok,stralf,dihstat,grmesh,ellcrd,anal_bs,cutoff_bs,residue_ref,early_vol,&
liganddef,ligand_id,readhet,nohet_ca,dihsave,eceppok,bsalign,motifpdb,motdl,ellpdb,plres_id,planetype,planedef,outpdb,final,&
savefit,overlay,bbaln,gridtype,kwsize,kwsize2,fold,aamax,estrtype,edihass,stranal,emesh,dccaprof,cmap,cmaptype,cmapcut,vear,&
saveseq,seqtype,num_mol,algn_typ,bcutoff,ecutoff,alnall,aadihsave,block,rtype,inpseqtype,AAoverlay,cmapmol,dorama,meanr,fit,&
dih_outliers,capri,recs,rece,ligs,lige)
write(16,*)'...parameters read in'

!00000000000000000000000000000000000000000000000000000000000000000000000
!some switches
if(rtype.eq.1)then
	write(16,*)'type of analysis: ',rtype,' - structural'
elseif(rtype.eq.2)then
	write(16,*)'type of analysis: ',rtype,' - sequential'
	goto 10
elseif(rtype.eq.3)then
	write(16,*)'type of analysis: ',rtype,' - structural and sequential'
else
	write(*,850)'rtype'
endif

!...controls for structure analysis     
CALL GETARG(2,BUFFER)
inquire(file=BUFFER,exist=fpresent)
if(.not.fpresent)then
	write(*,806)BUFFER
	stop
endif
READ(BUFFER,*)pdb1
write(16,*)'Input PDB file: ',pdb1

if(num_mol.gt.1)then
CALL GETARG(3,BUFFER)
inquire(file=BUFFER,exist=fpresent)
	if(.not.fpresent)then
		write(*,806)BUFFER
		stop
	endif
READ(BUFFER,*)pdb2
write(16,*)'Second input PDB file: ',pdb2
endif

!......controls for sequences analysis	
10	if(rtype.ne.1)then
  	 if(inpseqtype.eq.1)then
		write(16,*)'sequence alignment' 
		CALL GETARG(4,BUFFER)
		READ(BUFFER,*)inpalign
		inquire(file=inpalign,exist=fpresent)	
	 elseif(inpseqtype.eq.2)then
		write(16,*)'aminoacid sequences'
		CALL GETARG(4,BUFFER)
		READ(BUFFER,*)aasekw
		inquire(file=aasekw,exist=fpresent)
	 elseif(inpseqtype.eq.3)then
		write(16,*)'structural sequences'
		CALL GETARG(4,BUFFER)
		READ(BUFFER,*)assekw
		inquire(file=assekw,exist=fpresent)
	 elseif(inpseqtype.eq.4)then
		write(16,*)'aminoacid and structural sequences'
		CALL GETARG(4,BUFFER)
		READ(BUFFER,*)aasekw
		inquire(file=aasekw,exist=fpresent)
		CALL GETARG(5,BUFFER)
		READ(BUFFER,*)assekw
		inquire(file=assekw,exist=fpresent2)
		if(.not.fpresent2)then
        		write(*,806)'assekw'
        	 	stop
         	endif
	 else
		write(*,850)'inpseqtype'
		stop
	 endif
	 if(.not.fpresent)then
         	write(*,806)'sequence file'
         	stop
         endif	  	
       	endif		

!0000000000000000000000000000000000000000000000000000000000000000000000
!....analysis of sequences
if(rtype.ne.1)then
!---------------------------------------------------------------------	
!.....reads in the sequence alignment
	if(inpseqtype.eq.1)then
		write(16,*)'analysis of sequence alignment'
		open(7,file=inpalign,status='old',err=100)
		write(16,*)'reads in alignment...'
		call readalign(inpalign,kwsize,resid,wind2,nseqaln,nrestab,seq_motif2,algn_typ,nmot,bcutoff,&
		ecutoff,ignmtf)
!write(*,*)'mot',seq_motif(1,1,1),seq_motif(1,1,2),seq_motif(1,2,1),seq_motif(1,2,2)
		write(16,*)'...alignment read in'
		close(7)
!---------------------------------------------------------------------	 
!.....reads in sequence	 
	elseif(inpseqtype.gt.1)then
	  write(16,*)'analysis of AA and SA sequences'
	  write(16,*)'reads in sequences...'
	  call readseq(resnam_str,resid_str,nres_str,resnam_aa,resid_aa,nres_aa,inpseqtype,nasseq,aasekw,assekw)
	  write(16,*)'...sequences read in'
 	endif
endif	
if(rtype.eq.2)goto 45
	
!000000000000000000000000000000000000000000000000000000000000000000000
!....analysis of structures
write(16,*)'analysis of structures'
write(16,*)'number of molecules for mutual analysis: ',num_mol 
nseq=0
numf=0

20 infile=pdb1
lpoz2=index(pdb1,' ')
namcor=pdb1(1:(lpoz2-1))
if(num_mol.gt.1.AND.numf.eq.1)then
	infile=pdb2
	lpoz2=index(pdb2,' ')
	fnamcor=namcor
	namcor=pdb2(1:(lpoz2-1))
endif
inquire(file=infile,exist=fpresent)
if(.not.fpresent)then
	write(*,806)infile
       	stop
endif
lpoz5=index(infile,'/',back=.TRUE.)
if(lpoz5.ne.0)namcor=infile(lpoz5+1:(lpoz2-1))
kpos2=index(namcor,'.')
kpos=index(namcor,' ')
numf=numf+1 !numeracja wczytywanych plików

open(1,file=infile,status='old',err=100)
write(16,*)'molecule #',numf,'	',namcor 
if(num_mol.eq.1)then	!liczba bialek/lancuchow podawana jako parametr.
	if(nseq.eq.0)nseq=1
elseif(num_mol.gt.1)then
	if(nseq.eq.1)then
		nseq=2
	elseif(nseq.eq.0)then
		nseq=1
	elseif(nseq.ne.2)then
		write(*,*)'WARNING: files count incorrect'
	endif
else
	write(*,850)'num_mol'
	stop
endif

nat=0
nhet=0
nlig=0
nres=0
nch=1
ncach(0)=0	
natch(0)=0	
protein=.FALSE.
drop=.FALSE.
refine=.TRUE.
bezwodor=.FALSE.
het_present=.FALSE.
if(nohydr.eq.1)then
	write(16,*)'WARNING: NO hydrogen atoms will be taken, set nohydr=0 to take them'
	bezwodor=.TRUE.
elseif(nohydr.ne.0)then
	write(*,850)'nohydr'
	stop
endif

!----------------------------------------------------------------------
!.....reading structural data
write(16,*)'reads coordinates...'
call readstr(nat,nres,nch,natch,ncach,coord,calfc,hetco,nhet,nseq,protein,atomy,bocz,linia1,linia2,filetp,elipsa,bezwodor,resnam,&
fresnam,am,protch,nmr_ok,drop,nonhet,het_present,liganddef,ligand_id,ligand,nlig,readhet,nohet_ca,resid,plres_id,planecrd,&
planetype,planedef,nplane,numf,flin1,flin2,fatom,fbocz,fnres,fnat,fresid,chainid)
write(16,*)'...structure read in'

if(.NOT.refine)then
	write(16,*)'WARNING: no structure REFINEMENT'
	drop=.TRUE.
endif
if(.NOT.protein)then
	write(16,*)'WARNING: NOT protein molecule'
	drop=.TRUE.
endif
if(nres.lt.minres)then
	write(16,*)'WARNING: AA # less than: ',minres
	drop=.TRUE.
endif
if(drop)then
	write(16,*)'DROP WARNING: structure', namcor,' not analysed'
      	nseq=nseq-1
        drop=.FALSE.
	stop
endif

!...prepares temperatue factore replacement
if(inpseqtype.eq.1)then
	do i=1,numf
	 do k=1,nrestab(i)
	  do n=1,bocz(k)
	tfrepl(k,n)=wind2(i,k)
	  end do	 
	 end do
	end do
endif		
	
!...aa codes translation
A: do n=1,nres
!write(*,*)numf,n,resnam(n)
        do l=1,20
        if(resnam(n).eq.aminok(l))then
        aaseq(n)=AA(l)
!	write(*,*)numf,n,resnam(n),' ',aaseq(n)
	CYCLE A
        endif
	end do
end do A
if(numf.eq.1)then
B: do n=1,fnres
        do l=1,20
        if(fresnam(n).eq.aminok(l))then
        faaseq(n)=AA(l)
	CYCLE B
        endif
	end do
end do B
endif
!00000000000000000000000000000000000000000000000000000000000000000000000000	
!computational modules

!-------------------------------------------------------------------------	
!.....transformation of the reference system
if(trnsys.eq.1)then
	write(16,*)'transforms system...'
	 if(refsys.eq.1)then
	  if(nseq.eq.1)goto 40
	 endif
call transform(coord,calfc,hetco,nhet,lp,lpca,nseq,refsys,planecrd,nch,ncach,natch,atomy,bocz,nonhet)
	write(16,*)'...system transformed'
elseif(trnsys.ne.0)then
	write(*,850)'trnsys'
	stop
endif
!-------------------------------------------------------------------------	
!.......calculating fi, psi angles
40	if(torsje.eq.1)then
	write(16,*)'calculates backbone dihedrals...'
	rewind 1
call fikpsik(fipsi,minres,ncach,nch,atomy,resid,resnam,namcor,dihsave,dih_outliers)
	write(16,*)'...dihedrals calculated'
	elseif(torsje.ne.0)then
        write(*,850)'torsje'
        stop
	endif
!-------------------------------------------------------------------------
!.....dihedrals sorting	
if(dihstat.eq.1)then
	write(16,*)'sorts out dihedrals within given mesh...'
	if(torsje.ne.1)then
	        write(*,*)' '
        	write(*,*)'ERROR: dihedrals needed - set: "torsje" = 1 '
        	write(*,*)' '
        	stop
        endif
	call dihsort(fipsi,grmesh,ngrcv,nres,aadihsave)
	write(16,*)'...dihedrals sorted out'
elseif(dihstat.ne.0)then
        write(*,850)'dihstat'
        stop
endif
!------------------------------------------------------------------------
!.......elliptic analysis 
if(elipsa.eq.1)then
	write(16,*)'performs elliptic analysis...'
	if(torsje.ne.1)then
		write(*,*)' '
		write(*,*)'ERROR: dihedrals needed to perform elliptic analysis - set:"torsje" = 1 '
		write(*,*)' '
		stop
	endif
	call ellipse(fipsi,nres,elldih,resnam,tekrok,asseq,stralf,minres,ncach,nch,namcor,ellfipsi)
	write(16,*)'...elliptic analysis done'
elseif(elipsa.ne.0)then
        write(*,850)'elipsa'
        stop
endif
!-----------------------------------------------------------------------	
!.......structures creation
45 if(fold.eq.1)then
	write(16,*)'creates protein structure...'
!	write(*,*)(resnam_aa(1,k),k=1,4)
	call structure(atomy,bocz,ellatomy,resid,ncach,nch,minres,ellcrd,nres_str,edihass,estrtype,aamax,resnam_str,resnam_aa,&
	resid_aa,fipsi,ellfipsi,namcor)
	write(16,*)'...structure created'
elseif(fold.ne.0)then
        write(*,850)'fold'
        stop  
endif	
if(rtype.eq.2)goto 101	
!----------------------------------------------------------------------
!......comparison of fi,psi elipse and native maps
if(dorama.eq.1)then
	write(16,*)'compares fi,psi elipse and native maps...'
	call ramadist(fipsi,ellfipsi,nch,ncach,minres,namcor,resid,resnam,eceppok,emesh,kwsize2,tfrepl,bocz,nres,meanr)
	write(16,*)'...maps compared'
elseif(dorama.ne.0)then
      	write(*,850)'dorama'
      	stop
endif
!-----------------------------------------------------------------------
!......analysis of the box       
if(volfit.eq.1)then
	write(16,*)'Protein volume analysis, choosing type of box definition' 
!	if(numf.eq.1)open(4,file=(namcor(1:(kpos2-1))//'.vol'),status='unknown',err=100)
		if(volref.eq.0)then
		write(16,*)'volume of the box calculated using: CA ATOMS'  
		elseif(volref.eq.1)then
		write(16,*)'volume of the box calculated using: ALL ATOMS'  
		elseif(volref.eq.2)then
		write(16,*)'volume of the box calculated using: SIDE-CHAIN GEOMETRIC CENTER'  
		elseif(volref.eq.3)then
		write(16,*)'volume of the box calculated using: SIDE-CHAIN MASS CENTER'  
		elseif(volref.eq.4)then
		write(16,*)'volume of the box calculated using: RESIDUE GEOMETRIC CENTER'  
		elseif(volref.eq.5)then
		write(16,*)'volume of the box calculated using: RESIDUE MASS CENTER'  
		elseif(volref.ne.0)then
			write(*,850)'volref'
			write(16,850)'volref'
			stop
		endif   
      
	write(16,*)'calculating box volume...'
	call volume(namcor,early_vol,bocz,atomy,am,volref,nres,ellatomy)
	write(16,*)'...box volume calculated'
!	close(4)
elseif(volfit.ne.0)then
      	write(*,850)'volfit'
      	stop
endif
!-----------------------------------------------------------------------      
!.....binding site analysis
if(anal_bs.eq.1)then
	write(16,*)'analyses binding site of a molecule...'
	write(16,*)'WARNING: following parameters schould be selected: torsje=1' 
	write(16,*)'ellipse=1,stralf=1, liganddef=1, ligand_id'
	if(het_present)then
	write(16,*)'analysis with hetero atoms'
	open(12,file=(namcor(1:(kpos2-1))//'.sa.bs'),status='unknown')
	open(13,file=(namcor(1:(kpos2-1))//'.aa.bs'),status='unknown')	
	open(4,file=(namcor(1:(kpos2-1))//'.cnt'),status='unknown')
	endif
call analyse_bs(atomy,coord,calfc,nres,cutoff_bs,resid,residue_ref,nseq,am,bocz,asseq,het_present,&
ligand,nlig,resnam,numf,resmap,bscenter,motdl,namcor,bsalign,fnamcor,aaseq,savefit,ncach,nch,molbm)
     	if(het_present)then
	close(12)
	close(13)
	close(4)
	endif
	write(16,*)'...binding site analyzed'
elseif(anal_bs.ne.0)then
      	write(*,850)'anal_bs'
      	stop
endif       

!--------------------------------------------------------------------      
!.....creation of the grid
if(gridok.eq.1)then
	write(16,*)'performs grid analysis...'
	write(16,*)'WARNING: parameters required: anal_bs=1, pdbsave=1'
	if(nseq.eq.2) open(4,file=(namcor(1:(kpos-1))//'.grd'),status='unknown',err=100)
	call grid(coord,nres,nseq,kwezel,space,bscenter,ucenter,ngcat,kcenter,grdtab,atomy,bocz,resid,grdcon,&
	seed,agrid,gpk,gpn,gridtype,ligand,nlig)
     	if(nseq.eq.2)close(4)
     	write(16,*)'...grid analysis done'
elseif(gridok.ne.0)then
      	write(*,850)'gridok'
      	stop
endif
!--------------------------------------------------------------------
!.....analyses structures
if(stranal.eq.1)then
	write(16,*)'structural analysis...'
	call sanal(atomy,fatom,bocz,fbocz,nres,fnres,fnat,ncach,nch,nseq,namcor,fnamcor,dccaprof,cmap,cmaptype,cmapcut,overlay,&
	numf,bbaln,AAoverlay,nat,cmapmol,num_mol,resnam,resid,chainid,mostki,fit,capri,recs,rece,ligs,lige)
	write(16,*)'...analysis done'
elseif(stranal.ne.0)then
	write(*,850)'stranal'
     	stop
endif
!-------------------------------------------------------------------
!...calculates VeR
if(vear.eq.1)then
	write(16,*)'calculates V and R parameters...'
	do m=1,nch
 	 open(4,file=(namcor(1:(kpos-1))//'.'//sfx(m)//'.vr'),status='unknown')
	 i=0
	 ws=0
	 nr=ncach(m)-ncach(m-1)
 	 do k=ncach(m-1)+1,ncach(m)
       		i=i+1
      		ws(i,1)=atomy(k,1,2)
		ws(i,2)=atomy(k,1,3)
		ws(i,3)=atomy(k,1,4)
		i=i+1
        	ws(i,1)=atomy(k,5,2)
		ws(i,2)=atomy(k,5,3)
		ws(i,3)=atomy(k,5,4)
		i=i+1
        	ws(i,1)=atomy(k,2,2)
		ws(i,2)=atomy(k,2,3)
		ws(i,3)=atomy(k,2,4)
		i=i+1
		ws(i,1)=atomy(k,3,2)
		ws(i,2)=atomy(k,3,3)
		ws(i,3)=atomy(k,3,4)
		i=i+1
		ws(i,1)=atomy(k,4,2)
		ws(i,2)=atomy(k,4,3)
		ws(i,3)=atomy(k,4,4)
	 end do		
	call vr(ws,nr)
	close(4)
	end do
	write(16,*)'...VR calculated'
elseif(vear.ne.0)then
	write(*,850)'vear'
	stop
endif
!--------------------------------------------------------------------	
!...superimposes structures of multiple alignment blocks
if(block.eq.1)then
  if(num_mol.gt.1.AND.numf.gt.1)then
	if(inpseqtype.eq.1)then
	  write(16,*)'sequence alignment based fragments definition'
	  write(16,*)'motifs declared in input parameter file eill be overriden'
	elseif(inpseqtype.gt.1)then
	  write(16,*)'aligned protein fragments defined by user'
	  write(16,*)'motifdef=1 schould be declared'
	elseif(inpseqtype.ne.0)then
	  write(*,850)'inpseqtype'
	  stop
	endif
	  write(16,*)'superimposes and compares structures ...'
	  call blocker(fnamcor,namcor,atomy,fatom,ncach,nch,num_mol,seq_motif2,fnat,numf,resid,fresid,nmot,&
	  fnres,nres,alnall,aaseq,faaseq,ignmtf,AAoverlay)
	  write(16,*)'...structure superimposed'
  elseif(num_mol.gt.1.AND.numf.eq.1)then
	write(16,*)'first structure, waiting for the second one' 
  else
	write(*,850)'num_mol'
	stop
  endif
elseif(block.ne.0)then
  write(*,850)'block'
  stop
endif	
!--------------------------------------------------------------------
!.....writes pdb file
if(pdbsave.eq.1)then
	write(16,*)'writes structure in PDB format...'
	rewind 1
	call crd2pdb(coord,nseq,hetco,nonhet,ellatomy,bocz,linia1,linia2,nres,ellpdb,nch,ncach,namcor,motifpdb,resmap,&
	atomy,grdtab,gridok,outpdb,minres,tfrepl,fatom,estrtype,aasekw,resnam_aa,resnam_str,nres_str)
	write(16,*)'...structure written'
elseif(pdbsave.ne.0)then
      	write(*,850)'pdbsave'
      	stop
endif
!.........................................................................
!...writes sequences
if(saveseq.eq.1)then
	write(16,*)'saves sequences...'
	do n=1,nres
	 do l=1,20
        	if(resnam(n).eq.aminok(l))aaseq(n)=AA(l)
	 end do
	end do	
	if(seqtype.eq.2)then !....writes AA and SA sequence in columns
	 write(16,*)'writes AA and SA sequences in columns...'
	 open(14,file=(namcor(1:(kpos2-1))//'.seq'),status='unknown')	
	 DO n=1,nres
		write(14,818)resid(n),aaseq(n),asseq(n)
	 END DO
	elseif(seqtype.eq.1)then !....writes AA sequences in fasta format
	 write(16,*)'writes AA sequences in fasta format...'
	 open(14,file=(namcor(1:(kpos2-1))//'-AA.fasta'),status='unknown')
	 do j=1,nch
		write(14,821)">",namcor(1:(kpos2-1))," chain #: ",j,(aaseq(n), n=ncach(j-1)+1,ncach(j))
	 end do		 
	elseif(seqtype.eq.3)then !....writes AA and SA sequences in fasta format
	 write(16,*)'writes AA and SA sequences in fasta format...'
	 open(14,file=(namcor(1:(kpos2-1))//'-AA_SA.fasta'),status='unknown')
	 do j=1,nch
		write(14,821)">",namcor(1:(kpos2-1)),"AA chain #: ",j,(aaseq(n), n=ncach(j-1)+1,ncach(j))
		write(14,821)">",namcor(1:(kpos2-1)),"SA chain #: ",j,(asseq(n), n=ncach(j-1)+1,ncach(j))
	 end do
	elseif(seqtype.eq.4)then!....writes SA sequences in fasta format
	 write(16,*)'writes SA sequences in fasta format...'
	 open(14,file=(namcor(1:(kpos2-1))//'-SA.fasta'),status='unknown')
	 do j=1,nch
	  write(14,821)">",namcor(1:(kpos2-1)),"SA chain #: ",j,(asseq(n), n=ncach(j-1)+1,ncach(j))
	 end do
	else
		write(*,850)'seqtype'
		stop
	endif

close(14)
write(16,*)'...sequences saved'
elseif(saveseq.ne.0)then
	write(*,850)'saveseq'
	stop
endif

!END SECTION
!takes next protein from the list
write(16,*)'end of run for: ',namcor
write(16,*)'takes second file'
close(1)
if(numf.lt.num_mol)goto 20      !take next file
	
!00000000000000000000000000000000000000000000000000000000000000000
!final analysis

!-------------------------------------------------------------------	
!.....performs final analysis
if(final.eq.1)then
	write(16,*)'performs final analysis...'
	call finanal(ngrcv,grmesh,dihstat,grdcon,gridok,fnamcor,flin1,flin2,fatom,fbocz,fnres)
	write(16,*)'...final analysis done'
elseif(final.ne.0)then
      	write(*,850)'final'
      	stop
endif

!0000000000000000000000000000000000000000000000000000000000000000000
!...finishing section
write(16,*)'end of analysis'
goto 101
	  
806     format(' '/&
       '******************    ERROR!    ******************'/&
       ' '/&
        'File "',a,'" doesn''t exist' )
818	format(i7,1x,a1,1x,a1)
821	format(a1,a10,1x,a12,1x,i4,/,10000a1)

850     format('ERROR: Wrong value of parameter: ',a8,/&
       'restart with proper value set')

100   write(16,*)'ERROR opening unit'
101   call CPU_TIME(tstop)
write(16,*)'program run terminated normally'
time=tstop-tstart
write(16,'(A20,1X,f10.2)')'processor time used [s]:',time
      
end


