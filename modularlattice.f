C modularlattice.f
C May 2023 David S. Goodsell, provided as is for unrestricted use
C Reads coordinates of RNApol, ribosomes and SMC from generate_pol.f, and
C DNA/RNA chain lengths and connectivity from generate_connect.f
C Creates lattice model of DNA, nascent/mature mRNA, and nascent protein
C
C Largely hardwired for syn3A with simple RNApol, ribosomes, and SMC
C
C In tests with syn3A models, fails ~half of the time due to crowding!
C Try again with a different random seed.
C
C must be compiled with subroutines.f
C
C input parameters:
C syn3A.pdb              #input coordinate file with RNApol, ribosome, SMC templates
C syn3A.connect          #input file that defines DNA/RNA chain lengths and connections
C syn3A_lattice.pdb      #output lattice coordinates
C 0                      #supercoiling number of DNA beads per turn
C 1234                   #random seed
C 60.                    #cell radius (lattice units)
C
	integer*4 lattice(-100:100,-100:100,-100:100)
	integer*4 coords_out(3500,0:9000,3)
	real*4 rcoords_out(3500,0:5000,3)
	integer*4 icoords(0:90000,3),end(3),ends(2,3)
	integer*4 ecoords_in(0:90000,3),ecoords_out(0:90000,3)
	integer*4 segends(3500,2,3),segdna(3500,5,5)
	integer*4 nDNAsegments,DNAsegments_ends(3500,2)
	integer*4 DNAsegments_length(3500)
	integer*4 DNAsegments_nsupercoil_bp(3500)
	integer*4 idir_preferred(3500)
	real*4 pcoords_supercoil(0:90000,3)
	real*4 rcoords_full(0:90000,3)
	real*4 rvector(3)

	integer*4 nlength
	integer*4 nlength_connect(0:3500)
	real*4 rpersistence
	integer*4 branch(0:50,2)

	character*80 line
	character*3 name1,name2,chaintype(3)

	character*3 segname(3500)
	integer*4 segtype(3500),seglength(3500)
	integer*4 segchaintype(3500)
	integer*4 nsegindex(3500),segindex(3500,10)
	integer*4 seg_nsupercoil_bp(3500)
	integer*4 nsegchain(3500),segchain(3500,10,3)
	integer*4 seg_persistence(3500)
	real*4 rseg_persistence(3500)
	integer*4 segmoltype(3500,2)

	character*3 name(10)
	integer*4 type(10,4)
	integer*4 nendindex(10),endindex(10,4)
	integer*4 nchainindex(10),chainindex(10,4,2)		

	integer*4 nchain,chainend(3500)

	integer*4 molcoords(3500,100,3)
	integer*4 moltype(3500)
	character*30 molstuff(3500,10),segstuff(3500,10)
	character*2 atomname
	character segchainname(3500)
	integer*4 chainbreak(3500),nchainbreak

C define molecule and chain types
C currently hardwired for RIB,POL,SMC
C chain types: 1=DNA, 2=RNA, 3=protein

	ntype=4

	chaintype(1)="DNA"
	chaintype(2)="RNA"
	chaintype(3)="PRO"

C these blocks define the control points and parts of template chains that will be included in the model
	name(1)="POL"
	type(1,1)=1
	type(1,2)=1
	type(1,3)=2
	type(1,4)=0
	nendindex(1)=3
	endindex(1,1)=1
	endindex(1,2)=5
	endindex(1,3)=6
	endindex(1,4)=0
	nchainindex(1)=2
	chainindex(1,1,1)=1
	chainindex(1,1,2)=5
	chainindex(1,2,1)=6
	chainindex(1,2,2)=6

	name(2)="RIB"
	type(2,1)=2
	type(2,2)=2
	type(2,3)=3
	type(2,4)=0
	nendindex(2)=3
	endindex(2,1)=1
	endindex(2,2)=5
	endindex(2,3)=6
	endindex(2,4)=0
	nchainindex(2)=2
	chainindex(2,1,1)=1
	chainindex(2,1,2)=5
	chainindex(2,2,1)=6
	chainindex(2,2,2)=6

	name(3)="SMC"
	type(3,1)=1
	type(3,2)=1
	type(3,3)=0
	type(3,4)=0
	nendindex(3)=4
	endindex(3,1)=1
	endindex(3,2)=5
	endindex(3,3)=6
	endindex(3,4)=10
	nchainindex(3)=2
	chainindex(3,1,1)=1
	chainindex(3,1,2)=5
	chainindex(3,2,1)=6
	chainindex(3,2,2)=10

	name(4)="TER"
	type(4,1)=0
	type(4,2)=0
	type(4,3)=0
	type(4,4)=0
	nendindex(4)=0
	endindex(4,1)=0
	endindex(4,2)=0
	endindex(4,3)=0
	endindex(4,4)=0
	nchainindex(4)=0
	chainindex(4,1,1)=0
	chainindex(4,1,2)=0
	chainindex(4,2,1)=0
	chainindex(4,2,2)=0

	read(5,101) line
	open(1,file=line)
	write(6,*) "Input rib/pol positions: ",line
	read(5,101) line
	open(3,file=line)
	write(6,*) "Input connectivity file: ",line
	read(5,101) line
	open(2,file=line)
	write(6,*) "Output lattice nucleoid: ",line

	read(5,*) nsupercoil_bp_in
	write(6,*) "Supercoiling beads/supercoil: ",nsupercoil_bp_in
	read(5,*) nrand
	write(6,*) "Random seed: ",nrand
	read(5,*) rcell
	write(6,*) "Cell radius (lattice units):",rcell

	npersistence_DNA=10
	rpersistence_DNA=1.-1./float(npersistence_DNA)
	npersistence_RNA=1
	rpersistence_RNA=1.-1./float(npersistence_RNA)

	do i=1,nrand
	r=rand()
	enddo

C set up lattice with cell perimeter
	icellsq=int(rcell**2)
	do i=-100,100
	do j=-100,100
	do k=-100,100
	lattice(i,j,k)=0
	irsq=i*i+j*j+k*k
	if (irsq.gt.icellsq) lattice(i,j,k)=999
	enddo
	enddo
	enddo

C read coordinates for POL RIB SMC
	nmol=0
	icount=0
	do imol=1,3500
	moltype(imol)=0
	enddo

	write(2,201) "MODEL ",0

 100	read(1,101,end=99) line
	if (line(14:15).eq."CM") write(2,101) line 
	
C this finds a header line for each molecule
	if ((line(1:3).ne."END").and.
     &      (line(1:4).ne."ATOM")) then
	nmol=nmol+1
	do itype=1,ntype
	if (line(1:3).eq.name(itype)) moltype(nmol)=itype 
	enddo
	endif
C reset atomcount for next molecule
	if (line(1:3).eq."END") then
	 icount=0
	endif

	if (line(1:4).eq."ATOM") then
	 icount=icount+1
	  read(line,102) rx,ry,rz
	  ix=int(rx)
	  iy=int(ry)
	  iz=int(rz)

	  lattice(ix,iy,iz)=1
	  if (line(14:15).eq."CI") lattice(ix,iy,iz)=2

	  if ((line(14:16).ne."CM").and.
     &        (line(14:15).ne."CI")) then
	   molcoords(nmol,icount,1)=ix
	   molcoords(nmol,icount,2)=iy
	   molcoords(nmol,icount,3)=iz
	   read(line,110) molstuff(nmol,icount)
	  endif

C ATOM endif
	endif

	goto 100
 101	format(a80)
 102	format(30x,3f8.3)
 103	format(3x,8x,3i8)
 110	format(a30)
 99	continue

	write(2,201) "ENDMDL"

	write(6,*) "Read molecules ",nmol
	write(6,*)

C CONNECTIVITY INFORMATION

	nsegments=0
	nchain=1
	ntot=0
	nchainbreak=0

 130	read(3,101,end=139) line

	if (line(1:3).eq."DNA") then
	  ichaintype=1
	  goto 130
	endif
	if (line(1:3).eq."RNA") then
	  ichaintype=2
	  goto 130
	endif
	if (line(1:3).eq."PRO") then
	  ichaintype=3
	  goto 130
	endif
	if (line(1:3).eq."END") then
	  chainend(nchain)=nsegments
	  nchain=nchain+1
	  nchainbreak=nchainbreak+1
	  chainbreak(nchainbreak)=nsegments+1
	  goto 130
	endif

	if (line(1:3).ne."SEG") goto 130

	read(line,120) iseg,name1,mol1,index1,
     &                   name2,mol2,index2,length,
     &                  segchainname(nsegments+1)
	nsegments=nsegments+1
	if (iseg.ne.nsegments) iseg=nsegments

	do itype=1,ntype
	if (name1.eq.name(itype)) segmoltype(iseg,1)=itype
	if (name2.eq.name(itype)) segmoltype(iseg,2)=itype
	enddo

	segchaintype(iseg)=ichaintype
	segname(iseg)=chaintype(ichaintype)

C segmoltype=4 for TER

	if (segmoltype(iseg,1).ne.4) then

	id1=endindex(segmoltype(iseg,1),index1)
	ic1=(index1+1)/2
	ijcount=0
	do i=1,3
	 segends(iseg,1,i)=molcoords(mol1,id1,i)
	enddo

	ido1=chainindex(segmoltype(iseg,1),ic1,1)
	ido2=chainindex(segmoltype(iseg,1),ic1,2)
	ido1=min(ido1,10)
	ido2=min(ido2,10)
	do j=ido1,ido2
	  ijcount=ijcount+1
	  do i=1,3
	  segchain(iseg,ijcount,i)=molcoords(mol1,j,i)
	  enddo
	  segstuff(iseg,ijcount)=molstuff(mol1,j)
	enddo
	nsegchain(iseg)=ijcount

	endif

	if (segmoltype(iseg,2).ne.4) then

	id2=endindex(segmoltype(iseg,2),index2)
	ic2=(id2+1)/2
	ijcount=0
	do i=1,3
	 segends(iseg,2,i)=molcoords(mol2,id2,i)
	enddo

	endif

	seglength(iseg)=length
	
cwrite(6,*) "DNA segment input ",iseg,seglength(iseg),
c    &    (segends(iseg,1,j),j=1,3),(segends(iseg,2,j),j=1,3)


	if (ichaintype.eq.1) then
	   seg_nsupercoil_bp(iseg)=nsupercoil_bp_in
	   rseg_persistence(iseg)=rpersistence_DNA
	   seg_persistence(iseg)=npersistence_DNA
	else
	   seg_nsupercoil_bp(iseg)=0
	   rseg_persistence(iseg)=rpersistence_RNA
	   seg_persistence(iseg)=npersistence_RNA
	endif

	idir_preferred(iseg)=int(rand()*6.)


	goto 130

 120     format(3x,i8,2x,a3,i8,i4,2x,a4,i8,i4,i8,2x,a1)

 139	continue


C  CONNECTIVITY DEFINED
	nchainbreak=nchainbreak-1
	write(6,*) "Connectivity Defined"
	write(6,*) "Segments,Chains: ",nsegments,nchainbreak+1
	write(6,*)
	

cdo iseg=1,nsegments
cwrite(6,*) iseg,segname(iseg),(segends(iseg,1,i),i=1,3),
c    &      "|",(segends(iseg,2,i),i=1,3),"|",seglength(iseg)
cnc=nsegchain(iseg)
cdo ic=1,nc
cwrite(6,*) segstuff(iseg,ic),(segchain(iseg,ic,i),i=1,3)
cenddo
cwrite(6,*)
cenddo
cstop

	if (nsegments.gt.3500) then
	  write(6,*) "too many segments: ",nsegments
	  stop
	endif

C BIG LOOP OVER SEGMENTS
	do iseg=1,nsegments

C FIRST: add all segments that connect two points

	if ((segmoltype(iseg,1).ne.4).and.
     &      (segmoltype(iseg,2).ne.4)) then

	nseglength=seglength(iseg)
	nsupercoil_bp=DNAsegments_nsupercoil_bp(iseg)
	npersistence_in=seg_persistence(iseg)
	npersistence_in=1

	do i=1,3
	ends(1,i)=segends(iseg,1,i)
	ends(2,i)=segends(iseg,2,i)
	enddo


C TBD: pass preferred start direction using ierror variable
cierror=idir_preferred(iseg)
	ierror=int(rand()*6.)

C uninsulate ends
	do i=-1,1
	do j=-1,1
	do k=-1,1
	if (lattice(ends(1,1)+i,ends(1,2)+j,
     &     ends(1,3)+k).eq.2) then
           lattice(ends(1,1)+i,ends(1,2)+j,
     &     ends(1,3)+k)=0
	endif
	if (lattice(ends(2,1)+i,ends(2,2)+j,
     &     ends(2,3)+k).eq.2) then
           lattice(ends(2,1)+i,ends(2,2)+j,
     &     ends(2,3)+k)=0
	endif
	enddo
	enddo
	enddo

	do itry=1,10000
	nrand=int(rand()*5000)
	call connect(ends,lattice,npersistence_in,
     &     nseglength,nrand,nseglength_return,icoords,ierror)	

	nlength_connect(iseg)=nseglength_return
	if ((nseglength_return.gt.2).and.(ierror.eq.0).and.
     &      (nseglength_return.le.seglength(iseg))) goto 400

	 if (itry.eq.10000) then
	  write(6,*) "too many tries in connect",iseg
	  rd=sqrt(float((ends(1,1)-ends(2,1))**2+
     &            (ends(1,2)-ends(2,2))**2+
     &            (ends(1,3)-ends(2,3))**2))
	  write(6,*) "requested segment length ",
     &     seglength(iseg)," last try ",
     &     nseglength_return, " dist ",rd
	  stop
	 endif


	enddo
 400	continue

cwrite(6,*) chaintype(segchaintype(iseg)),
c    &  	" connection okay: ",iseg,nseglength_return,itry

	do iout=1,nseglength
	do i=1,3
	coords_out(iseg,iout,i)=icoords(iout,i)
	enddo
	enddo
	ntot=ntot+nseglength_return

C endif for type test for connecting segment
	endif
C iseg loop
	enddo


	write(6,*)
	write(6,*) "Finished two-molecule connections"
	write(6,*)

cgoto 888
C SECOND fill out full segment length for two-molecule connections

	do iseg=1,nsegments

	if ((segmoltype(iseg,1).ne.4).and.
     &      (segmoltype(iseg,2).ne.4)) then

	nseglength=seglength(iseg)
	nsupercoil_bp=seg_nsupercoil_bp(iseg)
	nlength_in=nlength_connect(iseg)
	npersistence_in=seg_persistence(iseg)

	if (nsupercoil_bp.eq.0) then
C ADDING EXTENSIONS

	do iatom=1,nlength_in
	do i=1,3
	ecoords_in(iatom,i)=coords_out(iseg,iatom,i)
	enddo
	enddo

	do itry=1,1000
	if (itry.eq.1000) then
	 write(6,*) "extension failed: ",iseg
	 stop
	endif
	call extend(ecoords_in,lattice,npersistence_in,
     &   nlength_in,nseglength,
     &    nrand,nlength_return,ecoords_out,ierror)

	if (ierror.eq.0) then
	do iout=1,nlength_return
	do i=1,3
	 rcoords_out(iseg,iout,i)=float(ecoords_out(iout,i))
	enddo
	enddo
	nlength_connect(iseg)=nlength_return
	ntot=ntot+nlength_return

cwrite(6,*) "added extensions: ",iseg,nlength_in,
c    &   nlength_return,itry

	goto 700
	endif

C itry loop
	enddo

	else


C ADDING PLECTONEME

C define plectoneme parameters
	nlength_plect=(nseglength-nlength_connect(iseg))/2+1
	iposition_root=nlength_connect(iseg)/2
	call define_plect(nlength_plect,nlength_root,ibranch,branch)

cwrite(6,*) "plec params ",nlength_plect, nlength_root,
c    &     ibranch,(branch(ib,1),ib=1,ibranch),
c    &     (branch(ib,2),ib=1,ibranch)

	do iatom=1,nlength_in
	do i=1,3
	rcoords_full(iatom,i)=float(coords_out(iseg,iatom,i))
	enddo
	enddo

	do i=1,3
	end(i)=coords_out(iseg,iposition_root,i)
	rvector(i)=coords_out(iseg,iposition_root+1,i)-
     &             coords_out(iseg,iposition_root-1,i)
	enddo
	rd=sqrt(rvector(1)**2+rvector(2)**2+rvector(3)**2)
	do i=1,3
	rvector(i)=rvector(i)/rd
	enddo

C  TRY TO ADD PLECTONEME
	do itry=1,20
	if (itry.eq.10) then
	iposition_root=iposition_root+1
	do i=1,3
	end(i)=coords_out(iseg,iposition_root,i)
	rvector(i)=coords_out(iseg,iposition_root+1,i)-
     &             coords_out(iseg,iposition_root-1,i)
	enddo
	rd=sqrt(rvector(1)**2+rvector(2)**2+rvector(3)**2)
	do i=1,3
	rvector(i)=rvector(i)/rd
	enddo
	endif


	if (itry.eq.20) then
	  write(6,*) "failed adding plectoneme: ",iseg 
	  stop
	endif

	nrand=int(rand()*5000)

cwrite(6,*) "starting plect subroutines ",iseg,nlength_plect,
c    &   nlength_root,(end(i),i=1,3),nrand

	call randomwalk(end,lattice,rpersistence_DNA,nlength_root,
     &      nrand,icoords,ierror)

	if (ierror.eq.0) then
cwrite(6,*) "finished randomwalk "
	else
cwrite(6,*) "random walk failed: ",itry
	goto 710
	endif

	call supercoil(icoords,nlength_root,nsupercoil_bp,rvector,
     &      pcoords_supercoil,ierror)
	call insert_plectoneme(rcoords_full,pcoords_supercoil,
     &      nlength_in,nlength_root,iposition_root,
     &      nlength_return,ierror)

	if (ibranch.gt.0) then
	do ib=1,ibranch

	iposition_branch=iposition_root+branch(ib,1)
	nlength_branch=branch(ib,2)

	do i=1,3
	end(i)=int(rcoords_full(iposition_branch,i))
	rvector(i)=rcoords_full(iposition_branch+1,i)-
     &             rcoords_full(iposition_branch-1,i)
	enddo
	rd=sqrt(rvector(1)**2+rvector(2)**2+rvector(3)**2)
	do i=1,3
	rvector(i)=rvector(i)/rd
	nlength_in=nlength_return
	enddo

cwrite(6,*) "starting branch ",iseg,ib,branch(ib,1),
c    &   branch(ib,2),(end(i),i=1,3),nrand

	call randomwalk(end,lattice,rpersistence_DNA,nlength_branch,
     &      nrand,icoords,ierror)
	call supercoil(icoords,nlength_branch,nsupercoil_bp,rvector,
     &      pcoords_supercoil,ierror)
	call insert_plectoneme(rcoords_full,pcoords_supercoil,
     &      nlength_in,nlength_branch,iposition_branch,
     &      nlength_return,ierror)

	enddo
	endif

	if (ierror.eq.0) then
	do iout=1,nlength_return
	do i=1,3
	 rcoords_out(iseg,iout,i)=rcoords_full(iout,i)
	enddo
	enddo
	nlength_connect(iseg)=nlength_return
	ntot=ntot+nlength_return
	goto 700
	endif

 710	continue
C itry loop
	enddo

C endif for nsupercoil_bp test
	endif

 700	continue



C test for connected segments
	endif
C big segment loop
	enddo

	write(6,*) "Finished filling out connected segments"
	write(6,*)
cstop

C THIRD add chains with free termini

	do iseg=1,nsegments

	if ((segmoltype(iseg,1).eq.4).or.
     &      (segmoltype(iseg,2).eq.4)) then

	nseglength_in=seglength(iseg)

	do i=1,3
	if (segmoltype(iseg,2).eq.4) then
	  end(i)=segends(iseg,1,i)
	else
	  end(i)=segends(iseg,2,i)
	endif
	enddo

	rpersistence=rseg_persistence(iseg)

cwrite(6,*) "end for random walk ",iseg,rpersistence,
c    &      nseglength_in,(end(i),i=1,3)

C remove insulation from mRNA end

	do i=-1,1
	do j=-1,1
	do k=-1,1
	if (lattice(end(1)+i,end(2)+j,end(3)+k).eq.2) then
            lattice(end(1)+i,end(2)+j,end(3)+k)=0
	endif
	enddo
	enddo
	enddo

	do itry=1,10000
	if (itry.eq.10000) then
	 write(6,*) "failed adding random walk ",iseg,nseglength_in
	 stop
	endif
	nrand=rand()*6000
	call randomwalk(end,lattice,rpersistence,nseglength_in,
     &     nrand,icoords,ierror)
	if (ierror.eq.0) then
c    write(6,*) "random walk okay: ",iseg,
c    &        seglength_in,itry
	    goto 701
	endif

cwrite(6,*) "random walk failed, trying again: ",ipol,itry

	enddo
 701	continue

	if (segmoltype(iseg,2).eq.4) then

	do iout=1,nseglength_in
	do i=1,3
	rcoords_out(iseg,iout,i)=float(icoords(iout,i))
	enddo
	enddo

	else

	do iout=1,nseglength_in
	do i=1,3
	rcoords_out(iseg,iout,i)=
     &     float(icoords(nseglength_in-iout+1,i))
	enddo
	enddo

	endif

C random walk test
	endif
C iseg loop
	enddo

cwrite(6,*) "random walk finished :",iseg,nlength_return

 50	format(a6,i8)
 200	format(a4,i7,2x,a9,i4,4x,3f8.3,2f6.2,4i8)

	write(6,*) "Finished unconnected segments"
	write(6,*)
	write(6,*) "Writing Coordinates"


C WRITE COORDINATES! 
	itot=0
	imodel=1
	isegout=0
	write(2,201) "MODEL ",1

	do iseg=1,nsegments
	itotseg=0

	do ib=1,nchainbreak
	  if (iseg.eq.chainbreak(ib)) then
	    imodel=imodel+1
	    isegout=0
	    itot=0
	    write(2,201) "ENDMDL"
	    write(2,201) "MODEL ",imodel
	  endif
	enddo

	do i=1,nsegchain(iseg)
	itot=itot+1
	itotseg=itotseg+1
	isegout=isegout+1
	write(2,300) "ATOM",itot,
     &    segstuff(iseg,itotseg)(14:15),
     &    chaintype(segchaintype(iseg)),
     &    segchainname(iseg),itotseg,
     &    (float(segchain(iseg,i,j)),j=1,3),
     &    2.,0.
	enddo
C temp factor = 0 for atoms that don't move

	atomname="P "
	if (segchaintype(iseg).eq.3) atomname="CA" 
	do i=1,seglength(iseg)
	itotseg=itotseg+1
	itot=itot+1
	isegout=isegout+1
	write(2,300) "ATOM",itot,
     &     atomname,
     &    chaintype(segchaintype(iseg)),
     &    segchainname(iseg),itotseg,
     &      (rcoords_out(iseg,i,j),j=1,3),1.,1.
	enddo
 300	format(a4,i7,2x,a2,2x,a3,a2,i4,4x,3f8.3,2f6.2,4i8)

 201	format(a6,i6)
	enddo
	write(2,201) "ENDMDL"

	end
