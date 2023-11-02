C generate_pol.f
C May 2023 David S. Goodsell, provided as is for unrestricted use
C reads ribosome positions and generates RNApol positions using a random walk
C SMC is placed between RNApol positions
C hardwired for lattice units of 10bp/bead
C
C must be compiled with subroutines.f
C
C input parameters:
C syn3Asmall_ribosomes.pdb     #ribosome positions (nm, should be centered on origin)
C 6c6u_mask.pdb                #mask file for RNApol
C 4v6k_mask.pdb                #mask file for ribosome
C 7nyx_mask.pdb                #mask file for SMC
C syn3A.pdb                    #output coordinate file (lattice units)
C 187                          #number of RNApol and SMC
C 201.26                       #radius of cell (nm)
C 1234                         #random seed
C 20.                          #max distance between successive RNApol
C RANDOM                       #RANDOM or CIRCLE placement of RNApol
C
	integer*4 map(-100:100,-100:100,-100:100)
	integer*4 rmap(-100:100,-100:100,-100:100)
	integer*4 polmask_in(500,3),ribmask_in(500,3)
	integer*4 polmask_out(500,3),ribmask_out(500,3)
	integer*4 smcmask_in(500,3),smcmask_out(500,3)
	integer*4 pol(5000,3),rib(5000,3),smc(5000,3)
	character*80 line
	character*30 polinfo(500),ribinfo(500),smcinfo(500)
	real*4 coords(0:10000,3)
		
C---------open files----------

C input file with ribosome positions in nm
	read(5,21) line
	open(3,file=line)
	write(6,*) "input file with ribosome positions: ",line

C input polymerase mask
	read(5,21) line
	open(10,file=line)
	write(6,*) "input mask file for polymerase: ",line
	npol_mask=1
 20	read(10,21,end=29) line
	if (line(1:4).eq."ATOM") then
	  read(line,22) polinfo(npol_mask),rx,ry,rz
	  polmask_in(npol_mask,1)=int(rx)
	  polmask_in(npol_mask,2)=int(ry)
	  polmask_in(npol_mask,3)=int(rz)
	  npol_mask=npol_mask+1
	endif
 21	format(a80)
 22	format(a30,3f8.3,2f6.2)
	goto 20
 29	continue
	npol_mask=npol_mask-1
	write(6,*) "read in polymerase mask: ",npol_mask

C input ribosome mask
	read(5,21) line
	open(11,file=line)
	write(6,*) "input mask file for ribosome: ",line
	nrib_mask=1
 30	read(11,21,end=39) line
	if (line(1:4).eq."ATOM") then
	  read(line,22) ribinfo(nrib_mask),rx,ry,rz
	  ribmask_in(nrib_mask,1)=int(rx)
	  ribmask_in(nrib_mask,2)=int(ry)
	  ribmask_in(nrib_mask,3)=int(rz)
	  nrib_mask=nrib_mask+1
	endif
	goto 30
 39	continue
	nrib_mask=nrib_mask-1
	write(6,*) "read in ribosome mask: ",nrib_mask

C input SMC mask
	read(5,21) line
	open(12,file=line)
	write(6,*) "input mask file for SMC: ",line
	nsmc_mask=1
 40	read(12,21,end=49) line
	if (line(1:4).eq."ATOM") then
	  read(line,22) smcinfo(nsmc_mask),rx,ry,rz
	  smcmask_in(nsmc_mask,1)=int(rx)
	  smcmask_in(nsmc_mask,2)=int(ry)
	  smcmask_in(nsmc_mask,3)=int(rz)
	  nsmc_mask=nsmc_mask+1
	endif
	goto 40
 49	continue
	nsmc_mask=nsmc_mask-1
	write(6,*) "read in SMC mask: ",nsmc_mask

C output file with full masks
	read(5,21) line
	open(1,file=line)
	write(6,*) "output file with full masks: ",line

C---------read parameters----------

	read(5,*) npol
	write(6,*) "number of polymerase: ",npol
	read(5,*) radius
	write(6,*) "cell radius (nm): ",radius
	read(5,*) irand
	write(6,*) "random seed: ",irand
	read(5,*) rmax
	write(6,*) "max distance between RNApol: ",rmax
	read(5,21) line
	icircle=0
	if (line(1:6).eq."CIRCLE") icircle=1
	write(6,*) "placement: ",line(1:6)
	write(6,*)

C	DNA
	
	do i=1,irand
	r=rand()
	enddo

C read ribosome positions
	nrib=1
 60	read(3,21,end=69) line
	if (line(1:4).eq."ATOM") then
	  read(line,62) rx,ry,rz
	   rib(nrib,1)=int(rx/3.4)
	   rib(nrib,2)=int(ry/3.4)
	   rib(nrib,3)=int(rz/3.4)
	  nrib=nrib+1
	endif
	goto 60
 62	format(30x,3f8.3)
 69	continue
	nrib=nrib-1
	write(6,*) "Ribosome positions input: ",nrib
	write(6,*)

C cell radius stuff
	radius=radius/3.4
	radiussq=radius**2
	write(6,*) "Cell radius (lattice units): ",radius
	write(6,*)
	rmaxsq=rmax**2

	do i=-100,100
	do j=-100,100
	do k=-100,100
	map(i,j,k)=0
	rmap(i,j,k)=0
	rsq=i*i+j*j+k*k
	if (rsq.gt.radiussq) then
	  map(i,j,k)=1
	  rmap(i,j,k)=1
	endif
	enddo
	enddo
	enddo

C PLACE RIBOSOMES FIRST IN KNOWN POSITIONS
	do irib=1,nrib

	do itry=1,10000
	if (itry.eq.10000) then
	  write(6,*) "Failed to place ribosome: ",irib
	  stop
	endif

	istate=int(rand()*18.)+1


C try to place at input position first, then try placing in increasingly large neighborhood
	if (itry.eq.1) then
	  ixjitter=0
	  iyjitter=0
	  izjitter=0
	endif
	if (itry.gt.1) then
	  ixjitter=int(rand()*5.)-2
	  iyjitter=int(rand()*5.)-2
	  izjitter=int(rand()*5.)-2
	endif
	if (itry.gt.100) then
	  ixjitter=int(rand()*15.)-7
	  iyjitter=int(rand()*15.)-7
	  izjitter=int(rand()*15.)-7
	endif
	if (itry.gt.1000) then
	  ixjitter=int(rand()*35.)-17
	  iyjitter=int(rand()*35.)-17
	  izjitter=int(rand()*35.)-17
	endif

	call rotate(ribmask_in,ribmask_out,nrib_mask,istate)

	do imask=1,nrib_mask
	  ix=ribmask_out(imask,1)+rib(irib,1)+ixjitter
	  iy=ribmask_out(imask,2)+rib(irib,2)+iyjitter
	  iz=ribmask_out(imask,3)+rib(irib,3)+izjitter
	  if (map(ix,iy,iz).ne.0) goto 70
	enddo

C	overlap test shows all okay at this point
Cwrite(6,*) "Placed ribosome: ",irib,itry
	write(1,71) "RIB",irib,(rib(irib,i),i=1,3),istate
	do imask=1,nrib_mask
	  ix=ribmask_out(imask,1)+rib(irib,1)+ixjitter
	  iy=ribmask_out(imask,2)+rib(irib,2)+iyjitter
	  iz=ribmask_out(imask,3)+rib(irib,3)+izjitter
	  write(1,73) ribinfo(imask)(1:22),irib,
     &     float(ix),float(iy),float(iz),1.,1.
	  map(ix,iy,iz)=2
	
	enddo
	  write(1,72) "ENDRIB"
	  goto 77
 71	format(a3,5i8)
 72	format(a6)
 73	format(a22,i4,4x,3f8.3,2f6.2)

 70	continue
C	itry loop
	enddo
 77	continue
C	irib loop
	enddo

C INITIAL PLACEMENT OF RNApol

	roffset=rand()*2.*3.14159
	ryangle=rand()*3.14159
	rzangle=rand()*3.14159*2.

	do itryconnect=1,1000
	do ipol=1,npol
	rangle=((2.*3.14159*float(ipol))/float(npol))+
     &     roffset
	do itry=1,100000
	if (itry.eq.100000) write(6,*)
     &        "initial pol placement failed: ",ipol

C----calculate trial positions for circle/random options
	if (icircle.eq.0) then
c----RANDOM placement
	rx=(rand()*(radius*2.+1.))-radius
	ry=(rand()*(radius*2.+1.))-radius
	rz=(rand()*(radius*2.+1.))-radius

	else
C----CIRCLE placement in an offset box that rotates once about an arbitrary axis
	rx1=(rand()*(radius*1.1+1.))-radius*0.1
	ry1=(rand()*(radius*0.5+1.))-radius*0.25
	rz1=(rand()*(radius*2.+1.))-radius

	rx2=rx1*cos(rangle)+ry1*sin(rangle)
	ry2=-rx1*sin(rangle)+ry1*cos(rangle)
	rz2=rz1

	rx3=rx2*cos(ryangle)+rz2*sin(ryangle)
	ry3=ry2
	rz3=-rx2*sin(ryangle)+rz2*cos(ryangle)

	rx=rx3*cos(rzangle)+ry3*sin(rzangle)
	ry=-rx3*sin(rzangle)+ry3*cos(rzangle)
	rz=rz3

	endif
C-----end RANDOM/CIRCLE option

	ispacing=6
	if (itry.gt.20000) ispacing=5
	if (itry.gt.40000) ispacing=4
	if (itry.gt.60000) ispacing=3
	if (itry.gt.80000) ispacing=0

	do i=-ispacing,ispacing
	do j=-ispacing,ispacing
	do k=-ispacing,ispacing
	itest=int(rx)+i
	jtest=int(ry)+j
	ktest=int(rz)+k
	if (map(itest,jtest,ktest).ne.0) goto 210
	enddo
	enddo
	enddo

	if (ipol.ne.1) then
	rtestsq=(rx-rxlast)**2+(ry-rylast)**2+(rz-rzlast)**2
	if (rtestsq.gt.rmaxsq) goto 210
	endif
	
	coords(ipol,1)=rx
	coords(ipol,2)=ry
	coords(ipol,3)=rz
	rxlast=rx
	rylast=ry
	rzlast=rz
 100	format(a4,i7,2x,a9,i4,4x,3f8.3,2f6.2)
	do i=-4,4
	do j=-4,4
	do k=-4,4
	itest=int(rx)+i
	jtest=int(ry)+j
	ktest=int(rz)+k
	if (map(itest,jtest,ktest).eq.0)
     &       map(itest,jtest,ktest)=9
	enddo
	enddo
	enddo

	goto 211

 210	continue

C itry loop
	enddo
 211    continue

C ipol loop
	enddo

	do i=-100,100
	do j=-100,100
	do k=-100,100
	  if (map(i,j,k).eq.9) map(i,j,k)=0
	enddo
	enddo
	enddo
	rconnect=sqrt((coords(1,1)-coords(npol,1))**2+
     &           (coords(1,2)-coords(npol,2))**2+
     &           (coords(1,3)-coords(npol,3))**2)
	write(6,*) "Pol connect distance:",itryconnect,rconnect
	if (rconnect.le.rmax) goto 299

C itryconnect loop
	enddo


 299	continue

C PLACE POLYMERASE MASKS
	do i=-100,100
	do j=-100,100
	do k=-100,100
	  if (map(i,j,k).eq.9) map(i,j,k)=0
	enddo
	enddo
	enddo

	do ipol=1,npol
	pol(ipol,1)=int(coords(ipol,1))
	pol(ipol,2)=int(coords(ipol,2))
	pol(ipol,3)=int(coords(ipol,3))

	do itry=1,1000
	if (itry.eq.1000) then
	   write(6,*) "Failed to place polymerase: ",ipol
	   stop
	endif

	istate=int(rand()*18.)+1

	if (itry.eq.1) then
	  ixjitter=0
	  iyjitter=0
	  izjitter=0
	endif
	if (itry.gt.1) then
	  ixjitter=int(rand()*5.)-2
	  iyjitter=int(rand()*5.)-2
	  izjitter=int(rand()*5.)-2
	endif
	if (itry.gt.100) then
	  ixjitter=int(rand()*15.)-7
	  iyjitter=int(rand()*15.)-7
	  izjitter=int(rand()*15.)-7
	endif

	call rotate(polmask_in,polmask_out,npol_mask,istate)

	do imask=1,npol_mask
	  ix=polmask_out(imask,1)+pol(ipol,1)+ixjitter
	  iy=polmask_out(imask,2)+pol(ipol,2)+iyjitter
	  iz=polmask_out(imask,3)+pol(ipol,3)+izjitter
	  if (map(ix,iy,iz).ne.0) goto 170
	enddo

C overlap test shows all okay at this point, write position
Cwrite(6,*) "Placed polymerase: ",ipol,itry
	write(1,71) "POL",ipol,(pol(ipol,i),i=1,3),istate
	do imask=1,npol_mask
	  ix=polmask_out(imask,1)+pol(ipol,1)+ixjitter
	  iy=polmask_out(imask,2)+pol(ipol,2)+iyjitter
	  iz=polmask_out(imask,3)+pol(ipol,3)+izjitter
	  write(1,73) polinfo(imask)(1:22),ipol,
     &     float(ix),float(iy),float(iz),1.,1.
	  map(ix,iy,iz)=1
	
	enddo
	  write(1,72) "ENDPOL"
	  goto 177

 170	continue
C	itry loop
	enddo
 177	continue
C	irib loop
	enddo

C PLACE SMC
C add one SMC between each pair of POL
	nsmc=npol
	do ismc=1,nsmc
	if (ismc.lt.nsmc) then
	smc(ismc,1)=(pol(ismc,1)+pol(ismc+1,1))/2
	smc(ismc,2)=(pol(ismc,2)+pol(ismc+1,2))/2
	smc(ismc,3)=(pol(ismc,3)+pol(ismc+1,3))/2
	else
	smc(ismc,1)=(pol(ismc,1)+pol(nsmc,1))/2
	smc(ismc,2)=(pol(ismc,2)+pol(nsmc,2))/2
	smc(ismc,3)=(pol(ismc,3)+pol(nsmc,3))/2
	endif

	do itry=1,1000
	if (itry.eq.1000) then
	   write(6,*) "Failed to place SMC: ",ismc
	   stop
	endif

	istate=int(rand()*18.)+1

	if (itry.eq.1) then
	  ixjitter=0
	  iyjitter=0
	  izjitter=0
	endif
	if (itry.gt.1) then
	  ixjitter=int(rand()*5.)-2
	  iyjitter=int(rand()*5.)-2
	  izjitter=int(rand()*5.)-2
	endif
	if (itry.gt.100) then
	  ixjitter=int(rand()*19.)-9
	  iyjitter=int(rand()*19.)-9
	  izjitter=int(rand()*19.)-9
	endif

	call rotate(smcmask_in,smcmask_out,nsmc_mask,istate)

	do imask=1,nsmc_mask
	  ix=smcmask_out(imask,1)+smc(ismc,1)+ixjitter
	  iy=smcmask_out(imask,2)+smc(ismc,2)+iyjitter
	  iz=smcmask_out(imask,3)+smc(ismc,3)+izjitter
	  if (map(ix,iy,iz).ne.0) goto 270
	enddo

C overlap test shows all okay at this point
	write(6,*) "Placed SMC: ",ismc,itry
	write(1,71) "SMC",ismc,(smc(ismc,i),i=1,3),istate
	do imask=1,nsmc_mask
	  ix=smcmask_out(imask,1)+smc(ismc,1)+ixjitter
	  iy=smcmask_out(imask,2)+smc(ismc,2)+iyjitter
	  iz=smcmask_out(imask,3)+smc(ismc,3)+izjitter
	  write(1,73) smcinfo(imask)(1:22),ismc,
     &     float(ix),float(iy),float(iz),1.,1.
	  map(ix,iy,iz)=1
	
	enddo
	  write(1,72) "ENDSMC"
	  goto 277

 270	continue
C	itry loop
	enddo

 277	continue
C	irib loop
	enddo
	
	end
