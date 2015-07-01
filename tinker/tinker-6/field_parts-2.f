c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2001 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c This program computes the electrostatic field at a point defined by the
c midpoint of two atoms.  This is written to specific usage in the Webb group
c for looking at the electrostatic field at the midpoint of a nitrile.
c Usage: ./field_parts <xyz> <CD atom index> <NE atom index> 
c <pair 1 starting atom index> <pair 1 ending atom index> <pair 2 starting
c atom index> <pair 2 ending atom index> ...
c
c
      program field_parts
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      include 'usage.i'
      include 'solute.i'
      include 'files.i'
      integer cdi,nei,ex0,ex1
      integer istart,iend
      real*8 bvmag
      real*8 extfield
      integer i,j,l,ia
      integer iz,ix
      real*8 a(3,3)
      logical exist,query
      character*120 string
      character*120 fieldfile
      character*7 ext
      integer ifield,iext,izeroes, freeunit
      real*8 cdpoint(3),nepoint(3)
      real*8 cdpi(8,3),nepi(8,3)
      real*8 nesigma(3), atom_site(100,3)
      integer excluded_atoms(100)
      real*8 potential(100)
      real*8 zvec(3),yvec(3),xvec(3)
      real*8 midpoint(3), bondvector(3)
      real*8 cdfield(3), nefield(3), tfield(3), exfield(3)
      real*8 field(3)
      real*8 cdf,nef,tf,exf
      real*8 ring_dist
      real*8 cdvecfield(3),nevecfield(3)
      real*8 cdsubtract(3),nesubtract(3)
      real*8 cx,cy,cz,nx,ny,nz
      integer nsites, nexcluded
      integer is_accounted_for
      nexcluded = 0
      nsites = 0
      ring_dist = 0.7d0
      nx = 0.0d0
      ny = 0.0d0
      nz = 0.0d0
      cx = 0.0d0
      cy = 0.0d0
      cz = 0.0d0
      do i = 1, 100
          excluded_atoms(i) = 0
          potential(i) = 0.0d0
          do j = 1, 3
             atom_site(i,j) = 0.0d0
          end do
      end do
      do i = 1, 3
          cdpoint(i) = 0.0d0
          nepoint(i) = 0.0d0
          zvec(i) = 0.0d0
          yvec(i) = 0.0d0
          xvec(i) = 0.0d0
          midpoint(i) = 0.0d0
          bondvector(i) = 0.0d0
          cdfield(i) = 0.0d0
          nefield(i) = 0.0d0
          tfield(i) = 0.0d0
          exfield(i) = 0.0d0
          field(i) = 0.0d0
          cdvecfield(i) = 0.0d0
          nevecfield(i) = 0.0d0
          cdsubtract(i) = 0.0d0
          nesubtract(i) = 0.0d0
      end do
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
      call nextarg (string,exist)
      write(iout,15) filename
   15 format('Reading ',a,/)                   
c
c    zero out then read the nitrile atom indices
c
      extfield = 0.0d0
      cdi = -1
      nei = -1
      ex0 = -1
      ex1 = -1
      if (exist) then
         read (string,*,err=10,end=10)  cdi
      else
         print *,'Usage: ./field_parts <CD Index> <NE Index> ',
     &                    '<Ignore contribution start> <Ignore ', 
     &                     'contribution end> <potential atom site 1> ',
     &                     '<potential atom site 2> <potential atom ',
     &                     'site 3> ... etc'
        print *,'NOTE: Ignore contributions and potential atom ',
     &                     'sites should NOT include CD and NE'
         stop
      endif
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nei
      else
         print *,'Usage: ./field_parts <CD Index> <NE Index> ',
     &                    '<Ignore contribution start> <Ignore ', 
     &                     'contribution end> <potential atom site 1> ',
     &                     '<potential atom site 2> <potential atom ',
     &                     'site 3> ... etc'
        print *,'NOTE: Ignore contributions and potential atom ',
     &                     'sites should NOT include CD and NE'
         stop
      endif
      call nextarg (string,exist)
      if (exist) then
          read (string,*,err=10,end=10) ex0
      endif
      call nextarg (string,exist)
      if (exist) then
          read (string,*,err=10,end=10) ex1
      endif
   10 continue
c
c     rotate the multipole components into the global frame
c
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
      end do
c
c     compute induced dipoles to be removed from QM multipoles
c
      call induce
c
c     if one of the PB implicit solvent model is used,
c     use the SCRF induced dipoles (uinds) instead of
c     vacuum induced dipoles (uind)
c
      if (solvtyp .eq. 'PB') then
         do i = 1,npole
            do l = 1, 3
                uind(l,i) = uinds(l,i)
            end do
         end do
      else if (solvtyp .eq. 'PB-HPMF') then
         do i = 1,npole
            do l = 1, 3
                uind(l,i) = uinds(l,i)
            end do
         end do
      end if
c
c     Find the CD-NE bond midpoint, bond vector, and magnitude of
c     the bond vector
c
      cdpoint(1) = x(cdi)
      cdpoint(2) = y(cdi)
      cdpoint(3) = z(cdi)
      nepoint(1) = x(nei)
      nepoint(2) = y(nei)
      nepoint(3) = z(nei)
      midpoint(1) = x(cdi)*.5d0 + x(nei)*.5d0
      midpoint(2) = y(cdi)*.5d0 + y(nei)*.5d0
      midpoint(3) = z(cdi)*.5d0 + z(nei)*.5d0
      bondvector(1) = x(nei) - x(cdi)
      bondvector(2) = y(nei) - y(cdi)
      bondvector(3) = z(nei) - z(cdi)
      call dotp(bondvector,bondvector,bvmag)
      bvmag = sqrt(bvmag)
c
c    Make the XYZ coordinate vectors about the nitrile
c
      do i = 1, 3
          zvec(i) = bondvector(i) / bvmag
          yvec(i) = 1.0d0
      end do
      if (zvec(1) .ne. 0.0d0) then
          yvec(1) = -(zvec(2)*yvec(2) + zvec(3)*yvec(3)) / zvec(1)
      else if (zvec(2) .ne. 0.0d0) then
          yvec(2) = -(zvec(1)*yvec(1) + zvec(3)*yvec(3)) / zvec(2)
      else if (zvec(3) .ne. 0.0d0) then
          yvec(3) = -(zvec(1)*yvec(1) + zvec(2)*yvec(2)) / zvec(3)
      end if
      call normalize(yvec)
      call crossp(yvec,zvec,xvec)
      call normalize(xvec)
      call ring_points(cdpoint,xvec,yvec,ring_dist,cdpi)
      call ring_points(nepoint,xvec,yvec,ring_dist,nepi)
      do i = 1, 3
          nesigma(i) = nepoint(i) + zvec(i) * ring_dist
          atom_site(1,i) = cdpoint(i)
          atom_site(2,i) = nepoint(i)
      end do
      excluded_atoms(1) = cdi
      excluded_atoms(2) = nei
      nsites = 2
      nexcluded = 2
c
c     Add induced dipole moments to the permanent dipoles
c
      do i = 1,npole
         rpole(2,i) = rpole(2,i) + uind(1,i)
         rpole(3,i) = rpole(3,i) + uind(2,i)
         rpole(4,i) = rpole(4,i) + uind(3,i)
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     Start writing the output to file
c
      ifield = freeunit()
      iext = 3
      izeroes = 0
      call numeral(izeroes,ext,iext)
      fieldfile = filename(1:leng)//'.'//ext(1:iext)//'E'
      call version (fieldfile,'new')
      open(unit=ifield,file=fieldfile,status='new')
c
c     CD field contribution
c
      call fieldpoint(midpoint,cdi,cdi,cdfield)
      call project(cdfield,bondvector,cdf)
      write(ifield,200,advance='no') ";",cdi
 200  format(a,i17)
c
c     NE field contribution
c
      call fieldpoint(midpoint,nei,nei,nefield)
      call project(nefield,bondvector,nef)
      write(ifield,300,advance='no') nei
 300  format(i18)

c
c     Total field
c
      call fieldpoint(midpoint,1,npole,tfield)
      call project(tfield,bondvector,tf)
      write(ifield,400,advance='no') 1," -",npole
 400  format(i10,a,i6)
      extfield = tf - cdf - nef
c
c     Remove the field contributions due to some atoms
c
      if ( ex1 .gt. -1 .and. ex0 .lt. npole ) then
          if ( ex1 .gt. npole ) then
              ex1 = npole
            end if
          call fieldpoint(midpoint,ex0,ex1,exfield)
          call project(exfield,bondvector,exf)
          write(ifield,500,advance='no') ex0," -",ex1
 500      format(i10,a,i6) 
          extfield = extfield - exf
          do i = ex0, ex1
              excluded_atoms(nexcluded+1) = i
              nexcluded = nexcluded + 1
          end do
      end if
c
c     Include the external field in the header
c
      write(ifield,550,advance='no') "External_Field"
 550  format(a18)
c
c     Add the CD and NE components to header
      write(ifield,555,advance='no') "X1","Y1","Z1","X2","Y2","Z2"
 555  format(a18,a18,a18,a18,a18,a18)
c
c     Get the field at the CD and NE atoms, less the excluded stuff
c
      call fieldpoint(cdpoint,1,npole,cdvecfield)
      call fieldpoint(nepoint,1,npole,nevecfield)
      do i = 1, nexcluded
          call fieldpoint(cdpoint,excluded_atoms(i),
     &                    excluded_atoms(i),cdsubtract)
          call fieldpoint(nepoint,excluded_atoms(i),
     &                    excluded_atoms(i),nesubtract)
          cdvecfield = cdvecfield - cdsubtract
          nevecfield = nevecfield - nesubtract
      end do
      call project(cdvecfield,zvec,cz)
      call project(cdvecfield,yvec,cy)
      call project(cdvecfield,xvec,cx)
      call project(nevecfield,zvec,nz)
      call project(nevecfield,yvec,ny)
      call project(nevecfield,xvec,nx)
c
c     Loop over the set of atoms indices for potential calculation
c
      do i = 1, npole
          istart = -1   
          call nextarg(string,exist)
          if (exist) then
              read (string,*,err=40,end=40) istart
          else
              exit
          endif
   40     continue
          atom_site(nsites+1,1) = x(istart)
          atom_site(nsites+1,2) = y(istart)
          atom_site(nsites+1,3) = z(istart)
          nsites = nsites + 1
c
c         Make sure the atom sites are in the excluded list
c         and add them if they are not
c
          do j = 1, nexcluded 
              if (istart .eq. excluded_atoms(j)) then 
                  is_accounted_for = 1
                  exit
              else
                  is_accounted_for = 0
              endif
              continue
         end do
         if (is_accounted_for .eq. 0) then
             excluded_atoms(nexcluded+1) = istart
             nexcluded = nexcluded + 1
         end if
      end do
      j = 1
c
c     Look over every atom site
c
      do i = 1, nsites
          write(ifield,560,advance='no') "p",j
 560      format(a17,i1)
          call potpoint(atom_site(i,:),excluded_atoms,nexcluded,
     &         potential(j))
          j = j + 1
      end do
c 
c     Look at the pi sites around CD
c
      write(ifield,570,advance='no') "pcpi",j
 570  format(a17,i1)
      do i = 1, 8
          call potpoint(cdpi(i,:),excluded_atoms,nexcluded,
     &         potential(j))
      end do
      j = j + 1
c
c     Look at the pi sites around NE
c
      write(ifield,580,advance='no') "pnpi",j
 580  format(a17,i1)
      do i = 1, 8
          call potpoint(nepi(i,:),excluded_atoms,nexcluded,
     &        potential(j))
      end do
      j = j + 1
c
c     Look at the sigma site along bond vector
c
      write(ifield,590,advance='no') "pnsg",j
 590  format(a17,i1)
      call potpoint(nesigma,excluded_atoms,nexcluded,
     &        potential(j))
      j = j + 1
c
c     Write the field values               
c
      write(ifield,600,advance='no') cdf,nef,tf
 600  format(/,3f18.8,3f18.8,3f18.8)
      if ( ex1 .gt. -1 .and. ex0 .lt. npole ) then
          write(ifield,700,advance='no') exf
 700      format(3f18.8)
      end if
      write(ifield,800,advance='no') extfield,cx,cy,cz,nx,ny,nz
 800  format(3f18.8,3f18.8,3f18.8,3f18.8,3f18.8,3f18.8,3f18.8)
c
c     Write the potential values
c
      do i = 1, j-1
          write(ifield,810,advance='no') potential(i)
 810      format(3f18.8)
      end do
      write(ifield,820,advance='no') 
 820  format(/)
c
c     rotate the multipoles from global back to local frame
c
      do i = 1,npole
         call rotmat (i,a)
         call invert (3,3,a)
         call rotsite (i,a)
      end do
c
c     copy the final multipoles back into local frame array
c
      do i = 1,npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine crossp --  compute the cross product of vectors ##
c     ##                                                             ##
c     #################################################################
c
      subroutine crossp( a, b, c )
      implicit none
      real*8 a(3),b(3),c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine normalize -- normalize a vector                 ##
c     ##                                                             ##
c     #################################################################
c
      subroutine normalize(a)    
      implicit none
      real*8 a(3)
      real*8 length,length2
      integer i
      length2 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
      length = sqrt(length2)
      do i = 1,3
          a(i) = a(i) / length
      end do 
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dotp -- dot product between 2 vectors size 3    ##
c     ##                                                             ##
c     #################################################################
c
      subroutine dotp(a, b, res)      
      implicit none
      real*8 a(3),b(3),res
      integer i
      res = 0.0d0
      do i = 1,3
          res = res + a(i)*b(i)
      end do 
      return 
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fieldpoint  --  electrostatic field at a point  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fieldpoint" calculates the electrostatic field at a grid
c     point "i" as the total electrostatic field of the monopole,
c     dipole + induced dipole, and quadrupole
c

      subroutine fieldpoint (ri, istart, iend, field )
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polpot.i'
      include 'polgrp.i'
      include 'couple.i'
      include 'polar.i'
      include 'units.i'
      include 'chgpot.i'
      real*8 ri(3), field(3)
      integer istart, iend, i, j, k
      real*8 r2,r,rr3,rr5,rr7
      real*8 pdotr
      real*8 qx,qy,qz,qdotr
      real*8 coulombs
      real*8 kbt
      real*8 permitivity_constant
      real*8 cfac
      real*8 mfld(3),dfld(3),qfld(3),rvec(3),qq(3)
c
c     To convert the fields from e/Angstrom**2 to KbT/eAngstrom
c
c     e0 is actually 1/(pi * permitivity_constant)
      permitivity_constant = 35950207149.4727056d0
c     coulombs/election
      coulombs = 6.2415096516d18
c     I'd like to use a user-input temperature, but am not sure
c     where that variable is in the Tinker source code?
      kbt = 298 * 1.380648813d-23
      cfac = 2.5d9 * permitivity_constant /
     &                   (coulombs * coulombs * kbt)
c     560.741404486965 versus 560.729000337724
      cfac =  electric / dielec / 0.5922d0
c
c     Set starting fields equal to zero
c
      do i = 1,3
          mfld(i) = 0.0d0
          dfld(i) = 0.0d0
          qfld(i) = 0.0d0
          rvec(i) = 0.0d0
          qq(i) = 0.0d0
          field(i) = 0.d0
      end do
c
c     Loop over the atoms we listed
c     Think about making the print out more readable
c
      if ( iend .gt. -1 .and. istart .lt. npole ) then
         if ( iend. gt. npole ) then
            iend = npole
         end if
c         print *,"Looking at atoms",istart,"to",iend
         do i = istart, iend
c
c     Distance information
c
            rvec(1) = ri(1) - x(i)
            rvec(2) = ri(2) - y(i)
            rvec(3) = ri(3) - z(i)
            call dotp(rvec,rvec,r2)
            r = sqrt(r2)
            if ( r2 .ne. 0 ) then
               rr3 = 1.0d0 / (r*r2)
               rr5 = 3.0d0 / (r*r2*r2)
               rr7 = 15.0d0 / (r*r2*r2*r2) 
            else
               rr3 = 0.0d0
               rr5 = 0.0d0
               rr7 = 0.0d0
            end if
c
c     Monopole part
c
            do j = 1, 3
                mfld(j) = mfld(j) + rr3 * pole(1,i) * rvec(j)  
            end do
c
c     Dipole part
c
            pdotr = pole(2,i)*rvec(1)+ pole(3,i)*rvec(2) + 
     &              pole(4,i)*rvec(3)
            do j = 1, 3
                dfld(j) = dfld(j) + rr5 * pdotr * rvec(j) - rr3 * 
     &                     pole(j+1,i)
            end do
c
c     Quadrupole part
c
c                qxx            qxy             qxz
            qq(1) = pole(5,i)*rvec(1) + 
     &              pole(6,i)*rvec(2) + 
     &              pole(7,i)*rvec(3)
c                qyx            qyy             qyz
            qq(2) = pole(6,i)*rvec(1) + 
     &              pole(9,i)*rvec(2) + 
     &              pole(10,i)*rvec(3)
c                qzx            qzy             qzz
            qq(3) = pole(7,i)*rvec(1) + 
     &              pole(10,i)*rvec(2) + 
     &              pole(13,i)*rvec(3)
            call dotp(qq,rvec,qdotr)
            do j = 1, 3
                qfld(j) = qfld(j) + rvec(j) * rr7 * qdotr - rr5 * qq(j) 
     &                    * 2.0d0
            end do
         end do
         do i = 1, 3
             field(i) = mfld(i) + dfld(i) + qfld(i)
             field(i) = field(i) * cfac
         end do
c
c     Project the field along the bond vector
c
c         call project(mfld,bvec,mf)
c         call project(dfld,bvec,df)
c         call project(qfld,bvec,qf)
c         mf = mf*cfac
c         df = df*cfac
c         qf = qf*cfac
c         print *,'BV :',bvx,bvy,bvz,bvmag,cfac
c         print *,'Quadrupole :',qfx,qfy,qfz
c         print *, mf, df, qf, mf+df+qf
c         print *," "
      end if
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine project -- project vector onto another vector   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine project( v1, v2, res )
      implicit none
      real*8 v1(3), v2(3), res
      real*8 v2length, v2length2, v1v2
      call dotp(v2,v2,v2length2)
      v2length = sqrt(v2length2)
      call dotp(v1,v2,v1v2)
      res = v1v2 / v2length
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ring_points -- find the 8 points in a ring      ##
c     ##                                                             ##
c     #################################################################
c
      subroutine ring_points( atomcenter, xvec, yvec, ring_dist, 
     &                        points )
      implicit none
      real*8 atomcenter(3), xvec(3), yvec(3), zvec(3)
      real*8 points(8,3)
      real*8 ring_dist, rsqrt2
      integer i,j
      rsqrt2 = 1.0d0 / sqrt(2.0d0)
      do i = 1, 3
          points(1,i) = atomcenter(i) + ring_dist * xvec(i)
          points(2,i) = atomcenter(i) - ring_dist * xvec(i)
          points(3,i) = atomcenter(i) + ring_dist * yvec(i)
          points(4,i) = atomcenter(i) - ring_dist * yvec(i)
          points(5,i) = atomcenter(i) + ring_dist * (xvec(i)+yvec(i)) 
     &                                            * rsqrt2
          points(6,i) = atomcenter(i) + ring_dist * (xvec(i)-yvec(i)) 
     &                                            * rsqrt2
          points(7,i) = atomcenter(i) - ring_dist * (xvec(i)+yvec(i)) 
     &                                            * rsqrt2
          points(8,i) = atomcenter(i) - ring_dist * (xvec(i)-yvec(i)) 
     &                                            * rsqrt2
      end do
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine potpoint -- electrostatic potential at a point  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "potpoint" calculates the electrostatic potential at a grid
c     point "i" as the total electrostatic field of the monopole,
c     dipole + induced dipole, and quadrupole
c

      subroutine potpoint(point, excluded_atom, nexcluded, potential) 
c
c These have been validated against the potpoint subroutine in potential.f!
c
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polpot.i'
      include 'polgrp.i'
      include 'couple.i'
      include 'polar.i'
      include 'units.i'
      include 'chgpot.i'
      real*8 point(3),potential
      integer nexcluded
      integer excluded_atom(100)
      integer i, j, k
      real*8 r2,r,rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 pdotr
      real*8 qx,qy,qz,qdotr
      real*8 coulombs
      real*8 kbt
      real*8 permitivity_constant
      real*8 cfac
      real*8 rvec(3),qq(3)
      real*8 monopole, dipole, quadrupole
      real*8 mfld(3),dfld(3),qfld(3),field(3)
      integer do_atom
c
c     To convert the fields from e/Angstrom**2 to KbT/eAngstrom
c
c     e0 is actually 1/(pi * permitivity_constant)
      permitivity_constant = 35950207149.4727056d0
c     coulombs/election
      coulombs = 6.2415096516d18
c     I'd like to use a user-input temperature, but am not sure
c     where that variable is in the Tinker source code?
      kbt = 298 * 1.380648813d-23
      cfac = 2.5d9 * permitivity_constant /
     &                   (coulombs * coulombs * kbt)
c     560.741404486965 versus 560.729000337724
      cfac =  electric / dielec / 0.5922d0
c
c     Set starting fields equal to zero
c
      do i = 1,3
          mfld(i) = 0.0d0
          dfld(i) = 0.0d0
          qfld(i) = 0.0d0
          rvec(i) = 0.0d0
          qq(i) = 0.0d0
          field(i) = 0.d0
      end do
      monopole = 0.0d0
      dipole = 0.0d0
      quadrupole = 0.0d0
c
c     Loop over the atoms
c
      do i = 1, npole
c
c     Distance information
c
         do_atom = 1
         do j = 1, nexcluded
             if ( i .eq. excluded_atom(j) ) then
                 do_atom = 0
                 exit
             end if
             continue
         end do
         if ( do_atom .eq. 1 ) then 
            rvec(1) = point(1) - x(i)
            rvec(2) = point(2) - y(i)
            rvec(3) = point(3) - z(i)
            call dotp(rvec,rvec,r2)
            r = sqrt(r2)
            if ( r2 .ne. 0 ) then
               rr1 = 1.0d0 / r
               rr3 = 1.0d0 / (r*r2)
               rr5 = 3.0d0 / (r*r2*r2)
            else
               print *,"distance zero",i
               rr1 = 0.0d0
               rr3 = 0.0d0
               rr5 = 0.0d0
            end if
c
c     Monopole part
c
            monopole = monopole + pole(1,i) * rr1
c
c     Dipole part
c
            pdotr = pole(2,i)*rvec(1)+ pole(3,i)*rvec(2) + 
     &              pole(4,i)*rvec(3)
            dipole = dipole + pdotr * rr3
c
c     Quadrupole part
c
c                qxx            qxy             qxz
            qq(1) = pole(5,i)*rvec(1) + 
     &              pole(6,i)*rvec(2) + 
     &              pole(7,i)*rvec(3)
c                qyx            qyy             qyz
            qq(2) = pole(6,i)*rvec(1) + 
     &              pole(9,i)*rvec(2) + 
     &              pole(10,i)*rvec(3)
c                qzx            qzy             qzz
            qq(3) = pole(7,i)*rvec(1) + 
     &              pole(10,i)*rvec(2) + 
     &              pole(13,i)*rvec(3)
            call dotp(qq,rvec,qdotr)
            quadrupole = quadrupole + qdotr*rr5
         end if
      end do
      monopole = monopole * cfac
      dipole = dipole * cfac
      quadrupole = quadrupole * cfac
      potential = potential + monopole + dipole + quadrupole
c      print *, potential, monopole, dipole, quadrupole
      return
      end
