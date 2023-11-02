C subroutines.f
C May 2022 David S. Goodsell, provided as is for unrestricted use
C-------------------------------------------------------------------
	subroutine connect(ends,lattice,npersistence,nlength_request,
     &    irand,nlength_return,coords_return,ierror)
C-------------------------------------------------------------------
C reads coordinates of two ends and a lattice of excluded points
C returns biased random walk that connects the endpoints

	integer*4 coords_return(0:90000,3),nlength_return
	integer*4 lattice(-100:100,-100:100,-100:100)
	integer*4 chain(0:90000,3),ends(2,3)
	integer*4 nlength,npersistence
	integer*4 idir(0:12,3)

c max number of tries when adding new structures
	nmaxstep=90000
	idirpreferred=ierror
	ierror=0
	nstep=0

	rpersistence2=1.-1./float(npersistence)

c set up random number generator
	do i=1,irand
	r=rand()
	enddo

	chain(0,1)=ends(1,1)
	chain(0,2)=ends(1,2)
	chain(0,3)=ends(1,3)
	lattice(ends(1,1),ends(1,2),ends(1,3))=1
	lattice(ends(2,1),ends(2,2),ends(2,3))=1

c defining directions of allowable steps
	do id1=0,12
	do id2=1,3
	idir(id1,id2)=0
	enddo
	enddo
	idir(0,3)=1
	idir(1,1)=-1
	idir(2,1)=1
	idir(3,2)=-1
	idir(4,2)=1
	idir(5,3)=-1
	idir(6,3)=1
	idir(7,1)=-1
	idir(8,1)=1
	idir(9,2)=-1
	idir(10,2)=1
	idir(11,3)=-1
	idir(12,3)=1

	idirlast=1
	idirtest=1
C pass in a preferred start direction using ierror variable
	if ((idirpreferred.ne.0).and.(idirpreferred.le.6)) then
	  idirlast=idirpreferred
	  idirtest=idirpreferred
	endif
	

C loop for adding points
	idmin=90000000
	idlast=9000000

	do istep=1,nmaxstep	

C fails because too many steps
	if (istep.eq.nmaxstep) then
	  ierror=1
	  do is=1,istep
	  lattice(chain(is,1),chain(is,2),chain(is,3))=0
	  enddo
	  return
	endif

	iran=int(rand()*6.)
	istepsuccess=0

	idmin=10000.
	idchainmax=0.
	icountmax=10

	do idir6=0,30

cif ((idir6.eq.0).and.(istep.ge.3)) then
	if (idir6.eq.0) then
C----persistence test, keep last direction?

           rpersistencetest=rand()
           if (rpersistencetest.lt.rpersistence2) then
             idirtest=idirlast
           else
             goto 200
           endif
	   ixs=chain(istep-1,1)+idir(idirtest,1)
	   iys=chain(istep-1,2)+idir(idirtest,2)
	   izs=chain(istep-1,3)+idir(idirtest,3)
	   idsq=((ixs-ends(2,1))**2+
     &           (iys-ends(2,2))**2+
     &           (izs-ends(2,3))**2)
	   if (idsq.eq.0) then
	     nstep=istep
             goto 100
	   endif

C----peristence test, else pick new direction
         else

	   iran=int(rand()*6.)
           idirtest=1+iran

	   ixs=chain(istep-1,1)+idir(idirtest,1)
	   iys=chain(istep-1,2)+idir(idirtest,2)
	   izs=chain(istep-1,3)+idir(idirtest,3)
	   idsq=((ixs-ends(2,1))**2+
     &          (iys-ends(2,2))**2+
     &           (izs-ends(2,3))**2)
	   if (idsq.eq.0) then
	     nstep=istep
             goto 100
	   endif

C----persistence test 
         endif

C check if the point is free
	   if (lattice(ixs,iys,izs).ne.0) goto 200
c check for dead ends
	     icount=0
	     do idirtest2=1,6
	      if (lattice(ixs+idir(idirtest2,1),iys+idir(idirtest2,2),
     &                    izs+idir(idirtest2,3)).ne.0) icount=icount+1
	     enddo
	     if (icount.ge.6) goto 200 
c if not a dead end, continue


	   if (idsq.lt.idlast) then
C if closer to end, automatically accept
	       ixmin=ixs
	       iymin=iys
	       izmin=izs
	       idmin=idsq
	       idlast=idsq
	       istepsuccess=1
	       idirlast=idirtest
	       goto 201

C else pick one with the emptiest neighborhood
	    else
	      if (icount.lt.icountmax) then
	       ixmin=ixs
	       iymin=iys
	       izmin=izs
	       idmin=idsq
	       idlast=idmin
	       idirlast=idirtest
	       icountmax=icount
               istepsuccess=2
	      endif
	   endif
C idir6 loop
 200	continue
	enddo

 201	continue

	if (istepsuccess.eq.0) then
	  ierror=1
	  do is=1,istep
	  lattice(chain(is,1),chain(is,2),chain(is,3))=0
	  enddo
	  return
	endif

cwrite(6,*) "final",ixmin,iymin,izmin,idirlast,idmin
cwrite(6,202) "ATOM",i,"PA  GLY C",i,
c    &      float(ixs),float(iys),float(izs),1.,1.

 202	format(a4,i7,2x,a9,i5,3x,3f8.3,2f6.2)
	chain(istep,1)=ixmin
	chain(istep,2)=iymin
	chain(istep,3)=izmin
	lattice(ixmin,iymin,izmin)=1

C	loop for adding points
	enddo

 100	continue
	nstep=nstep-1


 500	continue
C	OUTPUT COORDS

	nlength_return=nstep
	ilast=0
	do i=1,nstep
	coords_return(i,1)=chain(i,1)
	coords_return(i,2)=chain(i,2)
	coords_return(i,3)=chain(i,3)
	enddo
	return
	end	

C-------------------------------------------------------------------
	subroutine extend(coords_in,lattice,npersistence,
     &    nstep,nlength_request,
     &    irand,nlength_return,coords_return,ierror)
C-------------------------------------------------------------------
C reads coordinates of a chain and a lattice of excluded points
C returns a chain filled out with loop-like extensions to the requested length

	integer*4 coords_in(0:90000,3)
	integer*4 chain(0:90000,3),coords_return(0:90000,3)
	integer*4 nstep,nlength_request,nlength_return
	integer*4 lattice(-100:100,-100:100,-100:100)
	integer*4 nlength,npersistence
	integer*4 idir(0:12,3)

c max number of tries when adding new structures
	nmaxstep=90000
	ierror=0

cwrite(6,*) "starting extend1 ",nstep, nlength_request,irand,npersistence

	do i=1,nstep
	do j=1,3
	chain(i,j)=coords_in(i,j)
	enddo
	enddo

	rpersistence=float(npersistence)
	rpersistence2=0.90

c set up random number generator
	do i=1,irand
	r=rand()
	enddo

c defining directions of allowable steps
	do id1=0,12
	do id2=1,3
	idir(id1,id2)=0
	enddo
	enddo
	idir(0,3)=1
	idir(1,1)=-1
	idir(2,1)=1
	idir(3,2)=-1
	idir(4,2)=1
	idir(5,3)=-1
	idir(6,3)=1
	idir(7,1)=-1
	idir(8,1)=1
	idir(9,2)=-1
	idir(10,2)=1
	idir(11,3)=-1
	idir(12,3)=1

	idirlast=1
	idirtest=1

	do iextension=1,nmaxstep

	if (iextension.eq.nmaxstep) then
	  ierror=2
	  do il=1,istep
	  lattice(chain(il,1),chain(il,2),chain(il,3))=0
	  enddo
	  return
	endif

	rstepdirect=float(nstep)
C values of 2 and 4 insulate the ends
	istep=int(2.+rand()*(rstepdirect-4.))

	iran=int(rand()*6.)
	do idir6=1,6
	 idirtest=idir6+iran
	   ixs1=chain(istep,1)+idir(idirtest,1)
	   iys1=chain(istep,2)+idir(idirtest,2)
	   izs1=chain(istep,3)+idir(idirtest,3)
	   ixs2=chain(istep+1,1)+idir(idirtest,1)
	   iys2=chain(istep+1,2)+idir(idirtest,2)
	   izs2=chain(istep+1,3)+idir(idirtest,3)

	idist=abs(ixs1-ixs2)+abs(iys1-iys2)+abs(izs1-izs2)

	if (idist.gt.1) then
c  write(6,*) "idist error: ",istep, idist,nstep
c  write(6,601) (chain(istep,j),j=1,3),(chain(istepnext,k),k=1,3)
c  write(6,601) (i,i=1,nstep)
 601	format(60i4)
	endif

	if ((lattice(ixs1,iys1,izs1).eq.0).and.
     &      (lattice(ixs2,iys2,izs2).eq.0)) then
	
	
	ipersistence=npersistence+int((rand()-0.5)*rpersistence/4.)
	do i=1,ipersistence
	   if ((lattice(ixs1,iys1,izs1).eq.0).and.
     &         (lattice(ixs2,iys2,izs2).eq.0)) then

	     lattice(ixs1,iys1,izs1)=1
	     lattice(ixs2,iys2,izs2)=1
	     do imove=nstep,istep,-1
	      do j=1,3
	       chain(imove+2,j)=chain(imove,j)
	      enddo
	     enddo
	     chain(istep+1,1)=ixs1
	     chain(istep+1,2)=iys1
	     chain(istep+1,3)=izs1
	     chain(istep+2,1)=ixs2
	     chain(istep+2,2)=iys2
	     chain(istep+2,3)=izs2
	     istep=istep+1
	     nstep=nstep+2
	     if (nstep.ge.nlength_request) goto 500

	     ixs1=ixs1+idir(idirtest,1)
	     iys1=iys1+idir(idirtest,2)
	     izs1=izs1+idir(idirtest,3)
	     ixs2=ixs2+idir(idirtest,1)
	     iys2=iys2+idir(idirtest,2)
	     izs2=izs2+idir(idirtest,3)
	   else
	     goto 400
	   endif
C	persistence loop
	enddo

	goto 400
C	found a good direction
	endif

C	testing directions
	enddo	

C	adding loops
 400	continue
	enddo	

 500	continue
C	OUTPUT COORDS

	nlength_return=nstep
	ilast=0
	do i=1,nstep
	coords_return(i,1)=chain(i,1)
	coords_return(i,2)=chain(i,2)
	coords_return(i,3)=chain(i,3)
 200    format(a4,i7,2x,a9,i4,4x,3f8.3,2f6.2,4i8)
c       write(3,200) "ATOM",i,"C   POL A",i,
c    &      (float(coords_return(i,j)),j=1,3),1.,1.,0,0,0,0
	enddo
	return
	end	

C-------------------------------------------------------------------
	subroutine randomwalk(end,lattice,rpersistence,nlength_request,
     &    irand,coords_return,ierror)
C-------------------------------------------------------------------
C inputs coordinates of one end and a lattice of excluded points
C returns a random walk of points with specified persistence and length

	integer*4 coords_return(0:90000,3)
	integer*4 lattice(-100:100,-100:100,-100:100)
	integer*4 chain(0:90000,3),end(3)
	integer*4 nlength
	real*4 rpersistence
	integer*4 idir(0:12,3)

c max number of tries when adding new structures
	nmaxstep=90000
	ierror=0
	nstep=0

c set up random number generator
	do i=1,irand
	r=rand()
	enddo

	chain(0,1)=end(1)
	chain(0,2)=end(2)
	chain(0,3)=end(3)
	lattice(end(1),end(2),end(3))=1

c defining directions of allowable steps
	do id1=0,12
	do id2=1,3
	idir(id1,id2)=0
	enddo
	enddo
	idir(0,3)=-1
	idir(1,1)=-1
	idir(2,1)=1
	idir(3,2)=-1
	idir(4,2)=1
	idir(5,3)=-1
	idir(6,3)=1
	idir(7,1)=-1
	idir(8,1)=1
	idir(9,2)=-1
	idir(10,2)=1
	idir(11,3)=-1
	idir(12,3)=1

C loop for adding points
 
	ixs=1
	iys=0
	izs=0
	iderlast=1
	idertest=1

	do istep=1,nlength_request

	istepsuccess=0

	iran=int(rand()*6.)

	do idir6=0,6

	 if ((idir6.eq.0).and.(istep.ge.3)) then
	     rpersistencetest=rand()
	   if (rpersistencetest.lt.rpersistence) then
             idirtest=idirlast
	   else
	     goto 200
	   endif
	 else
	   idirtest=idir6+iran
	 endif

	   ixs=chain(istep-1,1)+idir(idirtest,1)
	   iys=chain(istep-1,2)+idir(idirtest,2)
	   izs=chain(istep-1,3)+idir(idirtest,3)
	   if (lattice(ixs,iys,izs).eq.0) then
c check for dead end
	     icount=0
	     do idirtest2=1,6
	      if (lattice(ixs+idir(idirtest2,1),iys+idir(idirtest2,2),
     &                    izs+idir(idirtest2,3)).ne.0) icount=icount+1
	     enddo
	     if (icount.lt.6) then 
	       chain(istep,1)=ixs
	       chain(istep,2)=iys
	       chain(istep,3)=izs
	       lattice(ixs,iys,izs)=istep
	       istepsuccess=1
	       idirlast=idirtest
	       goto 100
	     endif
	   endif
C   idir6 loop
 200	continue
	enddo

	if (istepsuccess.eq.0) then
cwrite(6,200) "ATOM",i,"PA  GLY C",i,
c    &      float(ixs),float(iys),float(izs),1.,1.
c200	format(a4,i7,2x,a9,i5,3x,3f8.3,2f6.2)
	  ierror=1
c  write(6,*) "random walk fails at step ",istep
	  do i=1,istep
	   lattice(chain(i,1),chain(i,2),chain(i,3))=0
	  enddo
	  return
	endif

C	loop for adding points
 100	continue
	enddo

C output step coordinates and lattice map update
	do i=1,nlength_request
	do j=1,3
	coords_return(i,j)=chain(i,j)
	enddo
	enddo

	return
	end	

C-------------------------------------------------------------------
	subroutine supercoil(coords,nlength,nsupercoil_bp,rvector,
     &      coords_return,ierror)
C-------------------------------------------------------------------
C inputs a chain of lattice coordinates and returns an off-lattice superhelix

	integer*4 coords(0:90000,3)
	real*4 coords_return(0:90000,3)
	real*4 pnext(0:90000,3),pnormal(0:90000,3),pv(0:90000,3)
	real*4 rvector(3)

	rsupercoil=360./float(nsupercoil_bp)
	rs=0.3

C assign dummy coordinates for 0 and nlength+1--needed for supercoiling vector calculation
	do j=1,3
	coords(nlength+1,j)=coords(nlength,j)*2-coords(nlength-1,j)
	coords(0,j)=coords(1,j)*2-coords(2,j)
	enddo

c -- offset for first step --
	do j=1,3
	  pnext(0,j)=float(coords(2,j)-coords(1,j))
	  pnormal(0,j)=float(coords(2,j)-coords(1,j))
	  pv(0,j)=0.
	enddo
	
C --sort out the cross-plectoneme vector for the first position in the supercoil
	if ((rvector(1).ne.0).or.(rvector(2).ne.0).or.
     &      (rvector(3).ne.0)) then
	  rv=sqrt(rvector(1)**2+rvector(2)**2+rvector(3)**2)
	  pv(0,1)=rvector(1)/rv
	  pv(0,2)=rvector(2)/rv
	  pv(0,3)=rvector(3)/rv
	  if ((pv(0,1).eq.(abs(pnext(0,1)))).and.
     &        (pv(0,2).eq.(abs(pnext(0,2)))).and.
     &        (pv(0,3).eq.(abs(pnext(0,3))))) then
c      write(6,*) "supercoiling vector parallel with next vector"
c      write(6,*) "permuting vector"
	      pv(0,1)=rvector(2)/rv
	      pv(0,2)=rvector(3)/rv
	      pv(0,3)=rvector(1)/rv
	  endif
	else
	  if (pnext(0,1).ne.0.) pv(0,2)=1.
	  if (pnext(0,2).ne.0.) pv(0,3)=1.
	  if (pnext(0,3).ne.0.) pv(0,1)=1.
	endif

c -- calculate offsets
	isharplast=0.

	do istep=1,nlength
	rl=0.
	rn=0.
	do j=1,3
	pnext(istep,j)=float(coords(istep+1,j)-coords(istep,j))
	pnormal(istep,j)=float(coords(istep+1,j)-coords(istep-1,j))
	rn=rn+(pnext(istep,j))**2
	rl=rl+(pnormal(istep,j))**2
	enddo
	do j=1,3
	pnext(istep,j)=pnext(istep,j)/sqrt(rn)
	pnormal(istep,j)=pnormal(istep,j)/sqrt(rl)
	enddo
	enddo

	do is=1,nlength

	rparallel=(pnext(is-1,1)*pnext(is,1))+
     &            (pnext(is-1,2)*pnext(is,2))+
     &            (pnext(is-1,3)*pnext(is,3))
	if (rparallel.gt.0.01) goto 600

c --  sharp turn

	rdiff=(pv(is-1,1)*pnormal(is,1))+
     &        (pv(is-1,2)*pnormal(is,2))+
     &        (pv(is-1,3)*pnormal(is,3))
	if ((rdiff.gt.0.99).or.(rdiff.lt.-0.99)) then
	rxc=-pnext(is-1,2)*pnormal(is,3)+
     &       pnext(is-1,3)*pnormal(is,2)
	ryc=-pnext(is-1,3)*pnormal(is,1)+
     &       pnext(is-1,1)*pnormal(is,3)
	rzc=-pnext(is-1,1)*pnormal(is,2)+
     &       pnext(is-1,2)*pnormal(is,1)
	rn=sqrt(rxc*rxc+ryc*ryc+rzc*rzc)
	if (rn.eq.0) then
C --vectors are exactly parallel
c  write(6,*) "sharp turn 0 -- setting arbitrary normal "
	  if (pnext(is-1,1).ne.0.) then
	    ryc=1.
	  else
	    rxc=1.
	  endif
	  rn=1.
	endif
	rlast=1.
	if (isharplast.ne.0) rlast=-1.
	pv(is,1)=rlast*rxc/rn
	pv(is,2)=rlast*ryc/rn
	pv(is,3)=rlast*rzc/rn
	isharplast=1
	goto 601
	endif

 600	continue
	isharplast=0
	rxc=-pv(is-1,2)*pnormal(is,3)+pv(is-1,3)*pnormal(is,2)
	ryc=-pv(is-1,3)*pnormal(is,1)+pv(is-1,1)*pnormal(is,3)
	rzc=-pv(is-1,1)*pnormal(is,2)+pv(is-1,2)*pnormal(is,1)
	rn=sqrt(rxc*rxc+ryc*ryc+rzc*rzc)
	if (rn.eq.0) then
	   rxc=0.
	   ryc=0.
	   rzc=0.
	   if (pnext(is-1,1).ne.0) ryc=1.
	   if (pnext(is-1,2).ne.0) rzc=1.
	   if (pnext(is-1,3).ne.0) rxc=1.
	rn=1.
	endif
	pv(is,1)=rxc/rn
	pv(is,2)=ryc/rn
	pv(is,3)=rzc/rn

 601	continue

C	calculate pv vectors (is) loop
	enddo

c -- write first strand of plectoneme
c -- vectors are now twisted by 90 deg every step, change to superhelical density
	do i=1,nlength
	ra0=float(mod(i,4))*90.
	rai=mod(float(i)*rsupercoil,360.)
	ra=(rai-ra0)*3.1415/180.
	rcos=cos(ra)
	rsin=sin(ra)
	rx=pnormal(i,1)
	ry=pnormal(i,2)
	rz=pnormal(i,3)
	pvx=pv(i,1)*(rcos+rx*rx*(1.-rcos))+
     &      pv(i,2)*(rx*ry*(1.-rcos)-rz*rsin)+
     &      pv(i,3)*(rx*rz*(1.-rcos)+ry*rsin)
	pvy=pv(i,2)*(rcos+ry*ry*(1.-rcos))+
     &      pv(i,3)*(ry*rz*(1.-rcos)-rx*rsin)+
     &      pv(i,1)*(ry*rx*(1.-rcos)+rz*rsin)
	pvz=pv(i,3)*(rcos+rz*rz*(1.-rcos))+
     &      pv(i,1)*(rz*rx*(1.-rcos)-ry*rsin)+
     &      pv(i,2)*(rz*ry*(1.-rcos)+rx*rsin)
	pv(i,1)=pvx	
	pv(i,2)=pvy	
	pv(i,3)=pvz	

	enddo

	do i=1,nlength
	rox=pv(i,1)*rs
	roy=pv(i,2)*rs
	roz=pv(i,3)*rs
	rx=(float(coords(i,1))-rox)
	ry=(float(coords(i,2))-roy)
	rz=(float(coords(i,3))-roz)
	coords_return(i,1)=rx
	coords_return(i,2)=ry
	coords_return(i,3)=rz
	enddo

	iout=nlength
	do i=nlength,1,-1
	iout=iout+1
	rox=pv(i,1)*rs
	roy=pv(i,2)*rs
	roz=pv(i,3)*rs
	rx=(float(coords(i,1))+rox)
	ry=(float(coords(i,2))+roy)
	rz=(float(coords(i,3))+roz)
	coords_return(iout,1)=rx
	coords_return(iout,2)=ry
	coords_return(iout,3)=rz
	enddo

	return
	end
C-------------------------------------------------------------------
	subroutine insert_plectoneme(coords,pcoords,
     &      nlength,nplength,nposition,
     &      nlength_return,ierror)
C-------------------------------------------------------------------
C inserts coordinates of a plectoneme into another coordinate set

	real*4 coords(0:90000,3)
	real*4 pcoords(0:90000,3)
	real*4 coords_return(0:90000,3)
	real*4 coords_tmp(0:90000,3)

	do i=1,nposition
	do j=1,3
	coords_tmp(i,j)=coords(i,j)
	enddo
	enddo

	insert=nposition
	
	do i=1,nplength*2
	do j=1,3
	coords_tmp(insert+i,j)=pcoords(i,j)
	enddo
	enddo

	insert=nplength*2

	do i=nposition+1,nlength
	do j=1,3
	coords_tmp(insert+i,j)=coords(i,j)
	enddo
	enddo

	nlength_return=nlength+nplength*2

	do i=1,nlength_return
	do j=1,3
	coords(i,j)=coords_tmp(i,j)
	enddo
	enddo

	return
	end

C -------------------------------------------------------
	subroutine rotate(icoord_in,icoord_out,
     &    ncoord,istate_in)
C -------------------------------------------------------
C returns coordinates for 1 of 24 possible 90 deg rotations
C istate in index for 24 possible rotations

	integer index(24,3),sign(24,3)
	integer*4 icoord_in(500,3),icoord_out(500,3)
	integer*4 ncoord,istate

C around Z (x y z) (-y x z) (-x -y z) (y -x z)
	data index(1,:)/1,2,3/
	data index(2,:)/2,1,3/
	data index(3,:)/1,2,3/
	data index(4,:)/2,1,3/
	data sign(1,:)/1,1,1/
	data sign(2,:)/-1,1,1/
	data sign(3,:)/-1,-1,1/
	data sign(4,:)/1,-1,1/

C around Z + 90 deg X (x -z y) (-y -z x) (-x -z -y) (y -z -x)
	data index(5,:)/1,3,2/
	data index(6,:)/2,3,1/
	data index(7,:)/1,3,2/
	data index(8,:)/2,3,1/
	data sign(5,:)/1,-1,1/
	data sign(6,:)/-1,-1,1/
	data sign(7,:)/-1,-1,-1/
	data sign(8,:)/1,-1,-1/

C around Z - 90 deg X (x z -y) (-y z -x) (-x z y) (y z x)
	data index(9,:)/1,3,2/
	data index(10,:)/2,3,1/
	data index(11,:)/1,3,2/
	data index(12,:)/2,3,1/
	data sign(9,:)/1,1,-1/
	data sign(10,:)/-1,1,-1/
	data sign(11,:)/-1,1,1/
	data sign(12,:)/1,1,1/

C around Z + 90 deg Y (-z y x) (-z x -y) (-z -y -x) (-z -x y)
	data index(13,:)/3,2,1/
	data index(14,:)/3,1,2/
	data index(15,:)/3,2,1/
	data index(16,:)/3,1,2/
	data sign(13,:)/-1,1,1/
	data sign(14,:)/-1,1,-1/
	data sign(15,:)/-1,-1,-1/
	data sign(16,:)/-1,-1,1/

C around Z - 90 deg Y (z y -x) (z x y) (z -y x) (z -x -y)
	data index(17,:)/3,2,1/
	data index(18,:)/3,1,2/
	data index(19,:)/3,2,1/
	data index(20,:)/3,1,2/
	data sign(17,:)/1,1,-1/
	data sign(18,:)/1,1,1/
	data sign(19,:)/1,-1,1/
	data sign(20,:)/1,-1,-1/

C around Z + 180 deg X (x -y -z) (-y -x -z) (-x y -z) (y x -z)
	data index(21,:)/1,2,3/
	data index(22,:)/2,1,3/
	data index(23,:)/1,2,3/
	data index(24,:)/2,1,3/
	data sign(21,:)/1,-1,-1/
	data sign(22,:)/-1,-1,-1/
	data sign(23,:)/-1,1,-1/
	data sign(24,:)/1,1,-1/

	if ((istate_in.gt.24).or.(istate_in.lt.1)) then
	istate=1
	else
	istate=istate_in
	endif

	do i=1,ncoord
	icoord_out(i,1)=sign(istate,1)*
     &          (icoord_in(i,index(istate,1)))
	icoord_out(i,2)=sign(istate,2)*
     &          (icoord_in(i,index(istate,2)))
	icoord_out(i,3)=sign(istate,3)*
     &          (icoord_in(i,index(istate,3)))
	enddo


	return
	end
C-------------------------------------------------------------------
	subroutine define_plect(nptot,nroot,ibranch,branch)
C-------------------------------------------------------------------
c input nptot--length of plectoneme superhelix (half of sequence length)
c returns nroot=length of root plectoneme
c ibranch=number of branches
c branch(i,1)=position, branch(i,2)=length

C Boles,White,Cozzarelli(1990)213,931 each superhelical segment ~1000bp
C length of plectoneme segment ~500bp, ~50beads at 10bp/bead

	integer*4 nroot,ibranch,nbranch,branch(0:50,2)
	integer*4 nptot 

c calculate the number of plectoneme segments, needs to be odd
	isegment=nptot/50
	if (mod(isegment,2).eq.0) isegment=isegment+1
	isegment=max(isegment,1)
	nsegment=nptot/isegment

c calculate length of root plectoneme
	isegroot=isegment/2+1
	ibranch=isegment-isegroot
	nbranch=nsegment
	nroot=nptot-(ibranch*nbranch)

C define branches
C tricky because need to define positions on a growing plectoneme

	branch(0,1)=0
	branch(0,2)=0

C odd branches on 5' side of plectoneme
	ibcount=1
	ninsert=nsegment
	if (ibranch.gt.0) then
	do ib=1,ibranch,2
	  branch(ibcount,2)=nbranch
	  branch(ibcount,1)=ninsert
	  ibcount=ibcount+1
	  ninsert=ninsert+nsegment*4
	enddo
	endif

C even branches on 3' side of plectoneme
	ninsert=ninsert+nsegment
	if (ibranch.gt.1) then
	do ib=2,ibranch,2
	  branch(ibcount,2)=nbranch
	  branch(ibcount,1)=ninsert
	  ibcount=ibcount+1
	  ninsert=ninsert+nsegment*4
	enddo
	endif

	return
	end

