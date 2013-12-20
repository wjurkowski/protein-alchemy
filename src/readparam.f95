SUBROUTINE readparam(filetp,nohydr,trnsys,refsys,nplane,torsje,elipsa,elldih,tekrok,volfit,volref,kcenter,gridok,ngcat,&
kwezel,space,ucenter,pdbsave,nonhet,mostki,minres,nmr_ok,stralf,dihstat,grmesh,ellcrd,anal_bs,cutoff_bs,&
residue_ref,early_vol,liganddef,ligand_id,readhet,nohet_ca,dihsave,eceppok,bsalign,motifpdb,motdl,ellpdb,plres_id,&
planetype,planedef,outpdb,final,savefit,overlay,bbaln,gridtype,kwsize,kwsize2,fold,aamax,estrtype,edihass,stranal,emesh,&
dccaprof,cmap,cmaptype,cmapcut,vear,saveseq,seqtype,num_mol,algn_typ,bcutoff,ecutoff,alnall,aadihsave,block,rtype,inpseqtype,&
AAoverlay,cmapmol,dorama,meanr,fit,dih_outliers,capri,recs,rece,ligs,lige)
USE tablice
IMPLICIT none

INTEGER::anal_bs,bsalign,cmap,cmaptype,dccaprof,dihsave,edihass,algn_typ,elipsa,elldih,ellcrd,ellpdb,eceppok,meanr,fit,&
early_vol,estrtype,filetp,final,fold,gridok,gridtype,kcenter,lmot,kwsize,kwsize2,liganddef,mostki,motdl,motifpdb,minres,ngcat,&
nmr_ok,nohydr,nonhet,nohet_ca,outpdb,motdef,overlay,bbaln,pdbsave,planedef,planetype,readhet,residue_ref,refsys,savefit,dihstat,&
stralf,stranal,trnsys,torsje,vear,volfit,volref,nplane(4,2),kwezel(3),kwezelx,kwezely,kwezelz,l,saveseq,seqtype,tekrok,num_mol,&
ecutoff(20),bcutoff(20),alnall,aadihsave,block,rtype,inpseqtype,AAoverlay,dorama,cmapmol,dih_outliers,capri,recs,rece,ligs,lige

REAL(8)::cutoff_bs,cmapcut,emesh,grmesh,r1,r2,spacex,spacey,spacez,ucenterx,ucentery,ucenterz,space(3),ucenter(3),aamax(20,14)
CHARACTER::ligand_id*3,plres_id*3

NAMELIST /basics/rtype,inpseqtype,algn_typ,filetp,minres,nmr_ok,nonhet,nohet_ca,nohydr,nplane,num_mol,planedef,planetype,&
plres_id,readhet,refsys,stranal,trnsys,anal_bs,gridok,dorama,elipsa,vear,torsje,final,volfit,volref,dihstat
NAMELIST /geomet/mostki,cmap,cmaptype,cmapcut,dccaprof,block,bcutoff,ecutoff,alnall,overlay,bbaln,AAoverlay,cmapmol,fit,&
capri,recs,rece,ligs,lige
NAMELIST /rama/meanr,kwsize2,eceppok,emesh
NAMELIST /ramastat/grmesh,aadihsave
NAMELIST /eliptic/early_vol,stralf,tekrok,estrtype,edihass,fold
NAMELIST /grid/gridtype,kcenter,kwezelx,kwezely,kwezelz,kwsize,ngcat,spacex,spacey,spacez,ucenterx,ucentery,ucenterz
NAMELIST /site/motdl,bsalign,cutoff_bs,liganddef,ligand_id,residue_ref
NAMELIST /inpout/dihsave,dih_outliers,ellpdb,elldih,ellcrd,motifpdb,outpdb,pdbsave,savefit,saveseq,seqtype

!default values for all parameters
!basics
stranal=0
anal_bs=0
gridok=0
elipsa=0
dorama=0
dihstat=0
algn_typ=1
filetp=1
inpseqtype=2
minres=50
nmr_ok=0
nohydr=0
nonhet=0
nohet_ca=0
nplane(1,1)=0
nplane(2,1)=0
nplane(3,1)=0
nplane(4,1)=0
nplane(1,2)=0
nplane(2,2)=0
nplane(3,2)=0
nplane(4,2)=0
num_mol=1
planedef=0
planetype=0
plres_id='RES'
readhet=0
refsys=2
rtype=1
trnsys=0
vear=0
torsje=0
final=0
volfit=0
volref=0
block=0

!geomet
AAoverlay=0
bbaln=0
overlay=0
fit=0
cmap=0
cmaptype=1
cmapmol=2
cmapcut=12.0
dccaprof=0
mostki=0
block=0
do l=1,20
bcutoff(l)=0
ecutoff(l)=0
end do
alnall=0
capri=0
recs=1
rece=1
ligs=1
lige=1

!rama
eceppok=0
emesh=5.0
meanr=0
kwsize2=5

!ramastat
grmesh=5.0
aadihsave=0

!eliptic
early_vol=0
stralf=0
tekrok=1
estrtype=1
edihass=0
fold=0

!grid
gridtype=2
kcenter=2
kwezelx=30
kwezely=30
kwezelz=30
kwsize=5
ngcat=1
spacex=5.0
spacey=5.0
spacez=5.0
ucenterx=0.0
ucentery=0.0
ucenterz=0.0
ucenter(1)=ucenterx
ucenter(2)=ucentery
ucenter(3)=ucenterz 
space(1)=spacex   
space(2)=spacey  
space(3)=spacez 
kwezel(1)=kwezelx 
kwezel(2)=kwezely 
kwezel(3)=kwezelz 

!site
motdl=7
bsalign=0
cutoff_bs=10.0
liganddef=0
ligand_id='UNK'
residue_ref=0

!inpout
dihsave=0
dih_outliers=0
elldih=0
ellcrd=0
ellpdb=0
motifpdb=0
outpdb=0
pdbsave=0
savefit=0
saveseq=0
seqtype=1

do l=1,14
	aamax(1,l)=alamax(l)
	aamax(2,l)=argmax(l)
	aamax(3,l)=asnmax(l)
	aamax(4,l)=aspmax(l)
	aamax(5,l)=cysmax(l)
	aamax(6,l)=glnmax(l)
	aamax(7,l)=glumax(l)
	aamax(8,l)=glymax(l)
	aamax(9,l)=hismax(l)
	aamax(10,l)=ilemax(l)
	aamax(11,l)=leumax(l)
	aamax(12,l)=lysmax(l)
	aamax(13,l)=metmax(l)
	aamax(14,l)=phemax(l)
	aamax(15,l)=promax(l)
	aamax(16,l)=sermax(l)
	aamax(17,l)=thrmax(l)
	aamax(18,l)=trpmax(l)
	aamax(19,l)=tyrmax(l)
	aamax(20,l)=valmax(l)
end do

read(2,NML=basics)
if(stranal.eq.1)read(2,NML=geomet)
if(dorama.eq.1)read(2,NML=rama)
if(dihstat.eq.1)read(2,NML=ramastat)
if(elipsa.eq.1)read(2,NML=eliptic)
if(gridok.eq.1)read(2,NML=grid)
if(anal_bs.eq.1)read(2,NML=site)
read(2,NML=inpout)

100 return

END SUBROUTINE
