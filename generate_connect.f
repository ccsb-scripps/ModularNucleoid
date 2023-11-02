C generate_connect.f
C May 2023 David S. Goodsell, provided as is for unrestricted use
C generates a ".connect" file used by ModularLattice
C this version is hardwired for a syn3A genome with 54300 beads, lattice units of 10 bp/bead
C reads genome information from an input file with RNApol positions and nascent mRNA length
C SMC loops are placed between RNApol positions greater than a user defined minimum DNA segment length
C mature mRNA is chosen randomly from nascent mRNA length, and one ribosome is placed at the center 
C
C input parameters:
C syn3A.genome      #genome file with RNApol positions and mRNA length (bp units)
C syn3A.connect     #output .connect file
C 503,187,187       #number of ribosomes, RNApol, and SMC
C 120               #minimum DNA length to add SMC (lattice units)

	integer*4 poloffset,smcoffset,riboffset
	integer*4 seglength
	integer*4 rnarib(1000)
	integer*4 position(200),length(200),segment(200)
	integer*4 direction(200),translated(200)
	integer*4 ipolindex(200,2)
	character*80 line
	character chainID

	read(5,1) line
	open(10,file=line)
	write(6,*) "input genome file",line

	read(5,1) line
	open(1,file=line)
	write(6,*) "output connect file",line
 1	format(a80)

	read(5,*) nrib,npol,nsmc
	write(6,*) "ribosomes, polymerases, SMC: ",
     &    nrib,npol,nsmc

	read(5,*) seglength
	write(6,*) "Minimum DNA segment length for SMC: ",seglength

C read in genomic information
C hardwired for 10 bp/bead
C hardwired for syn3A genome with 54300 beads

	do i=1,npol
	read(10,11) position(i),length(i),segment(i),
     &   direction(i),translated(i)
	position(i)=position(i)/10
	length(i)=length(i)/10
	segment(i)=segment(i)/10
	if (direction(i).eq.1) then
	 ipolindex(i,1)=1
	 ipolindex(i,2)=2
	else
	 ipolindex(i,1)=2
	 ipolindex(i,2)=1
	endif
	enddo
 11	format(3i8,2i3)

	ntot=0
	do i=1,npol
	ntot=ntot+segment(i)
	enddo
	segment(npol+1)=segment(1)+54300-ntot
	ntot2=0
	do i=2,npol+1
	ntot2=ntot2+segment(i)
	enddo
	write(6,*)
	write(6,*) "total DNA length: ",ntot,ntot2

	poloffset=nrib
	smcoffset=nrib+npol
	riboffset=0

	write(1,10) "DNA",1
 10	format(a3,i8)
	icount=0

C DNA
C-------SMC loops with only DNA---------------

	itot=0

	do ipol=1,npol
C----
	if (segment(ipol+1).gt.seglength+40) then
C -5 accounts for segments inside masks
	ipolsmc=seglength/2-5
	ismc=segment(ipol+1)-seglength-5
	itot=itot+ipolsmc*2+ismc

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "POL",poloffset+ipol,ipolindex(ipol,2),
     &      "SMC",smcoffset+ipol,1,
     &      ipolsmc,"A"

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "SMC",smcoffset+ipol,2,
     &      "SMC",smcoffset+ipol,3,
     &      ismc,"A"

	ipolnext=ipol+1
	if (ipol.eq.npol) ipolnext=1
	icount=icount+1
	write(1,20) "SEG",icount,
     &      "SMC",smcoffset+ipol,4,
     &      "POL",poloffset+ipolnext,ipolindex(ipolnext,1),
     &      ipolsmc,"A"

C----
	else
	itot=itot+segment(ipol+1)
	ipolsmc=segment(ipol+1)/2-5

	icount=icount+1

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "POL",poloffset+ipol,ipolindex(ipol,2),
     &      "SMC",smcoffset+ipol,1,
     &      ipolsmc,"A"

	ipolnext=ipol+1
	if (ipol.eq.npol) ipolnext=1
	icount=icount+1
	write(1,20) "SEG",icount,
     &      "SMC",smcoffset+ipol,2,
     &      "POL",poloffset+ipolnext,ipolindex(ipolnext,1),
     &      ipolsmc,"A"


C----
	endif

	enddo

	write(1,20) "SEG",icount,
     &      "POL",poloffset+1,ipolindex(1,1),
     &      "TER",0,0,0,"A"
	write(1,20) "END"

 20	format(a3,i8,2x,a3,i8,i4,3x,a3,i8,i4,i8,2x,a1)

C nascent mRNA

	ichaincount=1

	do ipol=1,npol
	ichaincount=ichaincount+1
	icount=0

	chainID="B"
	if (translated(ipol).eq.0) chainID="C"

	write(1,10) "RNA",ichaincount

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "TER",0,0,
     &      "POL",poloffset+ipol,3,
     &      length(ipol),chainID

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "POL",poloffset+ipol,3,
     &      "TER",0,0,0,chainID
	write(1,20) "END"

	enddo

C mature mRNA

	do i=1,1000
	r=rand()
	enddo

	do irib=1,nrib
	ichaincount=ichaincount+1
	icount=0

 600	irna=int(rand()*float(npol))+1
	rnarib(irib)=length(irna)
	if (translated(irna).eq.0) goto 600

cwrite(6,*) "mature RNA half length ",irib,rnarib(irib)

	write(1,10) "RNA",ichaincount

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "TER",0,0,
     &      "RIB",riboffset+irib,1,
     &      rnarib(irib),"D"

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "RIB",riboffset+irib,2,
     &      "TER",0,0,rnarib(irib),"D"
	write(1,20) "END"

	enddo

C nascent protein

	do irib=1,nrib
	ichaincount=ichaincount+1
	icount=0

	write(1,10) "PRO",ichaincount

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "TER",0,0,
     &      "RIB",riboffset+irib,3,
     &      rnarib(irib)/3,"E"

	icount=icount+1
	write(1,20) "SEG",icount,
     &      "RIB",riboffset+irib,3,
     &      "TER",0,0,0,"E"
	write(1,20) "END"

	enddo

	end
