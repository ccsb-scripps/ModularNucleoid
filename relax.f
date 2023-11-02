C relax.f
C May 2023 David S. Goodsell, provided as is for unrestricted use
C Relaxes a set of lattice coordinates using a simple force field
C Molecule positions and control points from generate_pol.f are held rigid
C
C input parameters:
C syn3A_lattice.pdb    #input lattice coordinates
C syn3A_relax.pdb      #output relaxed coordinates
C 1500                 #optimization steps

	integer*4 lattice(-100:100,-100:100,-100:100)
	real*4 rlattice(-100:100,-100:100,-100:100,3)
	integer*4 neighbors(0:175000,0:2000)
	integer*4 nhist(0:100)
	real*4 coordsin(0:175000,4),coords(0:175000,4)
	real*4 maskcoords(500,3)
	real*4 rcutoff(0:175000)
	character*30 stuff(0:175000)
	real*4 delta(0:175000,3)
	integer*4 bonds(0:175000,10),nbonds(0:175000)
	integer*4 constraints(0:175000,5),nconstraints(0:175000)
	real*4 rconstraintlength(0:175000,5)
	integer*4 bondhist(-500:500),consthist(-500:500)
	integer*4 nonbondhist(0:1000)
	character*80 filename
	character*98 line
	integer*4 chainend(1500)

	read(5,20) filename
	open(1,file=filename)
	write(6,*) "input coords: ",filename
	read(5,20) filename
	open(2,file=filename)
	write(6,*) "output coords: ",filename
 20	format(a80)
	read(5,*) nstep
	write(6,*) "nstep: ",nstep


	do i=0,175000
	neighbors(i,0)=0
	enddo
	do i=0,100
	nhist(i)=0
	enddo

	rbondlength=1.0
	rbondweight=0.05
	rconstraintweight=0.025
	rnonbondweight=0.005
	rweightlattice=0.1
	randweight=0.001
	ircellsq=60*60
	rcutoff_DNA=2.
	rcutoff_RNA=1.5
	
	do i=-100,100
        do j=-100,100
        do k=-100,100

        lattice(i,j,k)=0
        irsq=i*i+j*j+k*k
        if (irsq.gt.ircellsq) lattice(i,j,k)=999

C rlattice contains weight-scaled vectors pointing to free space
	if (lattice(i,j,k).ne.0) then
	  ri=float(i)
	  rj=float(j)
	  rk=float(k)
	  rd=sqrt(ri*ri+rj*rj+rk*rk)
	  rlattice(i,j,k,1)=-ri/rd*rweightlattice
	  rlattice(i,j,k,2)=-rj/rd*rweightlattice
	  rlattice(i,j,k,3)=-rk/rd*rweightlattice
	else
	  rlattice(i,j,k,1)=0.
	  rlattice(i,j,k,2)=0.
	  rlattice(i,j,k,3)=0.
	endif

        enddo
        enddo
        enddo

c read and write model 0 coordinates for POL RIB SMC
	npoint=0
	imask=0
	imollast=1
 110	read(1,201,end=99) line
	write(2,201) line
	imask=imask+1
	read(line,111,end=99) imol,(maskcoords(imask,i),i=1,3)
 111	format(22x,i4,4x,3f8.3)

	if (imol.ne.imollast) then

c find center of each mask
	rxcenter=0
	rycenter=0
	rzcenter=0
	do i=1,imask-1
	rxcenter=rxcenter+maskcoords(i,1)/float(imask-1)
	rycenter=rycenter+maskcoords(i,2)/float(imask-1)
	rzcenter=rzcenter+maskcoords(i,3)/float(imask-1)
	enddo

c fill out rlattice vectors using centers of each mask
	do i=1,imask-1
	do i2=-1,1
	do j2=-1,1
	do k2=-1,1
	ix=int(maskcoords(i,1))+i2
	iy=int(maskcoords(i,2))+j2
	iz=int(maskcoords(i,3))+k2
	rxv=maskcoords(i,1)-rxcenter
	ryv=maskcoords(i,2)-rycenter
	rzv=maskcoords(i,3)-rzcenter
	rnorm=sqrt(rxv*rxv+ryv*ryv+rzv*rzv)
	if (rnorm.eq.0) rnorm=1.
	rlattice(ix,iy,iz,1)=(rxv/rnorm)*rweightlattice
	rlattice(ix,iy,iz,2)=(ryv/rnorm)*rweightlattice
	rlattice(ix,iy,iz,3)=(rzv/rnorm)*rweightlattice
	lattice(ix,iy,iz)=1
	npoint=npoint+1
	enddo
	enddo
	enddo
	enddo

	maskcoords(1,1)=maskcoords(imask,1)
	maskcoords(1,2)=maskcoords(imask,2)
	maskcoords(1,3)=maskcoords(imask,3)
	imask=1
	imollast=imol

	endif

	if (line(1:4).eq."ENDM") goto 120
	goto 110

 120	continue

	do i=-80,80,2
	write(6,2) (int(abs(rlattice(i,j,0,1)/rweightlattice*9.)),
     &     j=-80,80,2)
	enddo
 2	format(81i1)
	write(6,*)

C-------READ in coordinates of chains
	natom=1
	idna=0
	nchain=0
 100	read(1,201,end=99) line
 201	format(a98)

	if (line(1:3).eq."END") then

	nchain=nchain+1
	chainend(nchain)=natom

	if (idna.eq.0) then
	 nbonds(natom-1)=0
	  nconstraints(natom)=0
	  nconstraints(natom-1)=0
	else
C******** hardwired for one DNA circle *********
	 nbonds(natom-1)=0
	  nconstraints(natom-6)=1
	  nconstraints(natom-5)=1
	  nconstraints(natom-4)=1
	  nconstraints(natom-3)=1
	  nconstraints(natom-2)=0
	  nconstraints(natom-1)=0
	  nconstraints(natom)=0
	 endif
	  goto 100
	endif

	if (line(1:4).ne."ATOM") goto 100

	read(line,200,end=99) stuff(natom),(coordsin(natom,i),i=1,3),
     &             rocc,rtemp
	coords(natom,1)=coordsin(natom,1)
	coords(natom,2)=coordsin(natom,2)
	coords(natom,3)=coordsin(natom,3)
	coords(natom,4)=rtemp
	bonds(natom,1)=natom+1
	nbonds(natom)=1

C atom is DNA
	if (stuff(natom)(18:20).eq."DNA") then
	  idna=1

	  constraints(natom,1)=natom+2
	  rconstraintlength(natom,1)=2.
	  nconstraints(natom)=1
	  if (natom.gt.1) constraints(natom-1,1)=natom-1+2
	  if (natom.gt.1) rconstraintlength(natom-1,1)=2.
	  if (natom.gt.2) constraints(natom-2,1)=natom-2+2
	  if (natom.gt.2) rconstraintlength(natom-2,1)=2.
	  nconstraints(natom)=1
	  rcutoff(natom)=rcutoff_DNA

	  if (rtemp.eq.1.) then
	    constraints(natom,2)=natom+6
	    rconstraintlength(natom,2)=6.
	    nconstraints(natom)=2
	    rcutoff(natom)=rcutoff_DNA
	  endif

C atom is not DNA
	else
	  idna=0
	  nconstraints(natom)=0
	  rcutoff(natom)=rcutoff_RNA
	endif
	natom=natom+1
	goto 100
 200	format(a30,3f8.3,2f6.2,4i8)
 99	continue
	natom=natom-1
	write(6,*) "number of coordinates ",natom

C build neighbors list
	do i=1,natom-1
	do j=i+1,natom

	if (j-i.lt.2) goto 300
	rx=coords(i,1)-coords(j,1)
	if (abs(rx).gt.7) goto 300
	ry=coords(i,2)-coords(j,2)
	if (abs(ry).gt.7) goto 300
	rz=coords(i,3)-coords(j,3)
	if (abs(rz).gt.7) goto 300
	rsq=rx*rx+ry*ry+rz*rz
	if (rsq.lt.49) then
	  neighbors(i,0)=neighbors(i,0)+1
	  neighbors(i,neighbors(i,0))=j
	endif
 300	continue

	enddo
	enddo

	do i=1,natom-1
	nhist(neighbors(i,0)/10)=nhist(neighbors(i,0)/10)+1
	nhist(0)=nhist(0)+neighbors(i,0)
	enddo
cwrite(6,*) "Neighbors List"
cdo i=1,40
cwrite(6,301) i,(neighbors(i,j),j=1,neighbors(i,0))
cenddo
 301	format(i8,40i5)
	write(6,*)
	write(6,*) "Neighbors Histogram"
	write(6,301) (nhist(i),i=0,40)	
	write(6,*)
	
	do ib=-10,10
	bondhist(ib)=0
	consthist(ib)=0
	enddo
	do ib=0,1000
	nonbondhist(ib)=0
	enddo

C BIG LOOP OVER STEPS
	do istep=1,nstep
	write(6,*) "STEP ",istep


C start periodic housekeeping
	if (mod(istep,50).eq.0) then

	write(6,*) "istep ",istep
	write(6,33) "bonds ",(bondhist(j),j=-10,10)
	write(6,33) "const ",(consthist(j),j=-10,10)
	write(6,33) "nonbd ",(min(999,nonbondhist(j)),j=0,20)
	write(6,*)
 33	format(a6,21i4)

	do ib=-10,10
	bondhist(ib)=0
	consthist(ib)=0
	enddo

	do ib=0,1000
	nonbondhist(ib)=0
	enddo

C build neighbors list
	do i=1,natom-1
	neighbors(i,0)=0
	do j=i+1,natom

	if (j-i.lt.3) goto 310
	rx=coords(i,1)-coords(j,1)
	if (abs(rx).gt.7) goto 310
	ry=coords(i,2)-coords(j,2)
	if (abs(ry).gt.7) goto 310
	rz=coords(i,3)-coords(j,3)
	if (abs(rz).gt.7) goto 310
	rsq=rx*rx+ry*ry+rz*rz
	if (rsq.lt.49) then
	  neighbors(i,0)=neighbors(i,0)+1
	  neighbors(i,neighbors(i,0))=j
	endif
 310	continue
	enddo
	enddo

	do i=0,40
	nhist(i)=0
	enddo

	do i=1,natom-1
	nhist(neighbors(i,0)/10)=nhist(neighbors(i,0)/10)+1
	nhist(0)=nhist(0)+neighbors(i,0)
	enddo

	write(6,*) "Neighbors List"
	do i=1,40
	enddo
	write(6,*)
	write(6,*) "Neighbors Histogram"
	write(6,301) (nhist(i),i=0,40)	
	write(6,*)

C end periodic housekeeping
	endif

	do iatom1=1,natom
	delta(iatom1,1)=0.
	delta(iatom1,2)=0.
	delta(iatom1,3)=0.
	enddo

C BIG LOOP OVER ATOMS
	do iatom1=1,natom

	 rvx=0.
	 rvy=0.
	 rvz=0.

	  do ibonds=1,nbonds(iatom1)
	    iatom2=bonds(iatom1,ibonds)

	    rvxb=coords(iatom2,1)-coords(iatom1,1)
	    rvyb=coords(iatom2,2)-coords(iatom1,2)
	    rvzb=coords(iatom2,3)-coords(iatom1,3)

	    rd=sqrt(rvxb*rvxb+rvyb*rvyb+rvzb*rvzb)
	    rdiff=rd-rbondlength
	    rweight=rdiff*rbondweight
	    rweight=min(rbondweight*2.,rweight)
	    rweight=max(-rbondweight*2.,rweight)

	    if (rd.eq.0.) then
	     rd=1.
	    endif

	    rvxb=rvxb/rd*rweight
	    rvyb=rvyb/rd*rweight
	    rvzb=rvzb/rd*rweight
	    delta(iatom1,1)=delta(iatom1,1)+rvxb
	    delta(iatom1,2)=delta(iatom1,2)+rvyb
	    delta(iatom1,3)=delta(iatom1,3)+rvzb
	    delta(iatom2,1)=delta(iatom2,1)-rvxb
	    delta(iatom2,2)=delta(iatom2,2)-rvyb
	    delta(iatom2,3)=delta(iatom2,3)-rvzb
	
	if (mod(istep,50).eq.1) then
	 ihist=int(rdiff*10.)
	if ((ihist.le.10).and.(ihist.ge.-10)) bondhist(ihist)=
     &                                       bondhist(ihist)+1
	endif
	  enddo

	  do iconstraints=1,nconstraints(iatom1)
	    iatom2=constraints(iatom1,iconstraints)

	    rvxb=coords(iatom2,1)-coords(iatom1,1)
	    rvyb=coords(iatom2,2)-coords(iatom1,2)
	    rvzb=coords(iatom2,3)-coords(iatom1,3)

	    rd=sqrt(rvxb*rvxb+rvyb*rvyb+rvzb*rvzb)
	    rdiff=rd-rconstraintlength(iatom1,iconstraints)
	    rweight=rdiff*rconstraintweight
	    rweight=min(rconstraintweight,rweight)
	    rweight=max(-rconstraintweight,rweight)

	    if (rd.eq.0.) then
	     rd=1.
	    endif

	      rvxb=rvxb/rd*rweight
	      rvyb=rvyb/rd*rweight
	      rvzb=rvzb/rd*rweight
	      delta(iatom1,1)=delta(iatom1,1)+rvxb
	      delta(iatom1,2)=delta(iatom1,2)+rvyb
	      delta(iatom1,3)=delta(iatom1,3)+rvzb
	      delta(iatom2,1)=delta(iatom2,1)-rvxb
	      delta(iatom2,2)=delta(iatom2,2)-rvyb
	      delta(iatom2,3)=delta(iatom2,3)-rvzb

	if (mod(istep,50).eq.1) then
	 ihist=int(rdiff*10.)
	if ((ihist.le.10).and.(ihist.ge.-10)) consthist(ihist)=
     &                                       consthist(ihist)+1
	endif
	  enddo
	  

C Nonbonded contacts
	do inonbond=1,neighbors(iatom1,0)
	iatom2=neighbors(iatom1,inonbond)

	    rvxb=coords(iatom2,1)-coords(iatom1,1)
	    rvyb=coords(iatom2,2)-coords(iatom1,2)
	    rvzb=coords(iatom2,3)-coords(iatom1,3)

	    rd=sqrt(rvxb*rvxb+rvyb*rvyb+rvzb*rvzb)

	    if (rd.lt.rcutoff(iatom1)) then
	    if (rd.eq.0.) then
	     rd=1.
	    endif
	    rvxb=-rvxb/rd*rnonbondweight
	    rvyb=-rvyb/rd*rnonbondweight
	    rvzb=-rvzb/rd*rnonbondweight
	    delta(iatom1,1)=delta(iatom1,1)+rvxb
	    delta(iatom1,2)=delta(iatom1,2)+rvyb
	    delta(iatom1,3)=delta(iatom1,3)+rvzb
	    delta(iatom2,1)=delta(iatom2,1)-rvxb
	    delta(iatom2,2)=delta(iatom2,2)-rvyb
	    delta(iatom2,3)=delta(iatom2,3)-rvzb
	 endif
	if (mod(istep,50).eq.1) then
	ihist=int(rd*5.)
	if ((ihist.ge.0).and.(ihist.le.1000)) nonbondhist(ihist)=
     $        nonbondhist(ihist)+1
	endif

 600	  continue

C	iatom2 loop
	enddo

C add random motion
	rvx=rvx+(rand()-0.5)*randweight
	rvy=rvy+(rand()-0.5)*randweight
	rvz=rvz+(rand()-0.5)*randweight

	ix=int(coords(iatom1,1))
	iy=int(coords(iatom1,2))
	iz=int(coords(iatom1,3))

	if (lattice(ix,iy,iz).ne.0) then
	delta(iatom1,1)=delta(iatom1,1)+rlattice(ix,iy,iz,1)
	delta(iatom1,2)=delta(iatom1,2)+rlattice(ix,iy,iz,2)
	delta(iatom1,3)=delta(iatom1,3)+rlattice(ix,iy,iz,3)
	endif
	
C	BIG iatom1 LOOP
	enddo

	do iatom1=1,natom

	if (coords(iatom1,4).ne.0.) then

	coords(iatom1,1)=coords(iatom1,1)+delta(iatom1,1)
	coords(iatom1,2)=coords(iatom1,2)+delta(iatom1,2)
	coords(iatom1,3)=coords(iatom1,3)+delta(iatom1,3)

	endif
	enddo

C	BIG istep LOOP
	enddo
	imod=1

	write(2,901) "MODEL ",imod

        do i=1,natom
	do itest=1,nchain
	  if (chainend(itest).eq.i) then
	    write(2,901) "ENDMDL"
	 imod=imod+1
	 write(2,901) "MODEL ",imod
	endif
	enddo
	
	if (stuff(i)(14:15).eq."P5") stuff(i)(14:15)="P "
	if (stuff(i)(14:15).eq."P3") stuff(i)(14:15)="P "
	if (stuff(i)(14:15).eq."CN") stuff(i)(14:15)="CA"
	if (stuff(i)(14:15).eq."CC") stuff(i)(14:15)="CA"
        write(2,900) stuff(i),
     &      (coords(i,j),j=1,3),1.,coords(i,4)
 900    format(a30,3f8.3,2f6.2,4i8)
        enddo
	write(2,901) "ENDMDL"
 901	format(a6,i6)

	end
