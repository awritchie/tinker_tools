

      program calc_field
      use sizes
      use atoms
      use files
      use iounit
      use inform
      use mpole
      use units
      use polar
      use solute
      use limits
      use bound
      use polpot
      use potent
      use rigid 
      use titles
      use atomid
      integer i,j,ixyz,iind,ifld,frame,freeunit,trimtext
      integer lext, uindsize, uindn, next
      integer v1i, v2i
      integer exi(maxatm), nexatoms, exatom
      logical exist, douind, quit, opened
      character*1 letter
      character*7 ext
      character*120 record,string,xyzfile,fldfile,indfile
      character*120 tmpc
      integer uindtag(maxatm)
      character*120 uindname(maxatm)
      real*8 field(9), exfield(9), ifield(9)
      real*8 bondvector(3),prfield(3),prexfield(3)
      real*8 pfield, pexfield, v1field, v2field
      real*8 kfac
      real*8 cx,cy,cz,bondlength
      kfac = debye * 299.8391483043805 / 2.5852
c     
c     set up the structure and mechanics calculation
c     
      call initial
      call getxyz
      call mechanic
      write(iout,10) filename
  10  format('Reading ',a,/)
      if (.not.use_polar .or. poltyp.eq.'DIRECT') then
         print *, "not using mutual polarization"
         call fatal
      end if
c     
c     get addition input parameters
c     
      call nextarg (string,exist)
      if (.not. exist) then
         write(iout,20)
  20     format(/,' Calculate induced dipoles [Y]',
     &           /,' Use induced dipole file [filename]')
  30     continue
         write(iout,40)
  40     format(/,' Enter the desired induced dipole method',
     &              ' [Y, filename] :  ',$)
         read(input,50,err=30) string
  50     format(a120) 
      end if
c
c     set control flags
c
      douind = .false.
      if (string .eq. 'Y') douind = .true.
      if (string .eq. 'y') douind = .true.
      if (douind) then 
         write(iout,60)
  60     format(' Will calculate induced dipole moments')
      else 
         next = 1
         i = trimtext(string)
         indfile = string
         write(iout,70) indfile(1:i)
  70     format(' Will read induced dipole moments from ',a,".")
         inquire(file=indfile,exist=exist)
         if (.not. exist) then
            write (iout, 75)
  75        format(/,' Unable to find Induced Dipole File')
            call fatal
         end if
      end if
c
c     get atoms defining bond vector 
c
      call nextarg(string,exist)
      if (exist) then
         read(string,*,err=80) v1i
      end if
      if (.not. exist) then
  80     continue
         write(iout,90)
  90     format(' Index of atom the bond vector points AWAY from:  ')
         read(input,*,err=80) v1i
      end if 
      call nextarg(string,exist)
      if (exist) then
         read(string,*,err=100) v2i
      end if 
      if (.not. exist) then
  100    continue
         write(iout,110)
  110    format(' Index of atom the bond vector points TOWARDS:  ')
         read(input,*,err=100) v2i
      end if 
      write(iout,120) v1i, v2i
  120 format(' Bond vector points from atoms: ',i0,' to ',i0,'.')
c
c     get the indices of atoms not to include in the field
c 
      do i = 1, 100 
          exi(i) = 0
      end do
      exi(1) = v1i
      exi(2) = v2i
      nexatoms = 2
      do while (string .ne. 'E') 
         call nextarg(string,exist)
         call upcase(string)
         if (exist) then
            if (string .eq. 'E') then
               goto 150
            end if
            read(string,*,err=130) exatom
            if (.not. ANY(exi(1:nexatoms) == exatom)) then
               exi(nexatoms+1) = exatom
               nexatoms = nexatoms + 1
            end if
         end if
         if (.not. exist) then
  130       continue
            write(iout,140)
  140       format(' Index of atoms to exclude field contribution:  ',
     &              /,' (Enter [E] to end adding excluded atoms) ')
            read(input,*,err=130) string       
            call upcase(string)
            if (string .eq. 'E') then
               goto 150
            end if
            read(string,*,err=130) exatom
            if (.not. ANY(exi(1:nexatoms) == exatom)) then
               exi(nexatoms+1) = exatom
               nexatoms = nexatoms + 1
            end if
         end if 
      end do
  150 continue    
      if (nexatoms .gt. 0 ) then
         write(iout,160)
  160    format(' Excluding atoms :')
         do i = 1, nexatoms
            write(iout,170,advance='no') exi(i)
  170       format(2x,i0)
         end do
         write(iout,180)
  180    format(/)
      end if
c
c     reopen the coorinate file and read the first structure
c
      frame = 0
      ixyz = freeunit()
      xyzfile = filename
      call suffix(xyz,'xyz','old')
      open(unit=ixyz,file=xyzfile,status='old')
      rewind(unit=ixyz)
      call readxyz(ixyz)
c
c     Open input induced dipole file
c     
      if (.not.douind) then
         iind = freeunit()
         open(unit=iind,file=indfile,status='old')
      end if
c
c     Open output field file
c
      ifld = freeunit()
      fldfile = xyzfile(1:leng)//'.'//ext(1:iext)//'fld'
      call version(fldfile,'new')
      open(unit=ifld,file=fldfile,status='new')
c
c     Include a header
c
      write(ifld,185,advance='no') 'Units KbT/eA','frame',
     &                'uind','total'
  185 format('#',a15,/,'#',a8,2x,a16,2x,a16,2x)
      do i = 1, nexatoms
         write(ifld,186,advance='no'), 'monopole-',exi(i),
     &                                 'dipole-',exi(i),
     &                                 'quadrupole-',exi(i)
  186    format(a15,i0,2x,a15,i0,2x,a15,i0,2x)
      end do
      write(ifld,187), 'less_perm_excl', 'less_excl'
  187 format(a16,2x,a16)
c      write(iout,186) 'Units KbT/eA','frame',
c     &                'uind','less_excluded','total'
c  186 format('#',a15,/,'#',a8,2x,a16,2x,a16,2x,a16)
c
c     perform analysis for each successive coodinate structure
c
      do while(.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            if (douind) then
               write (iout,190) 
  190          format(' Exting after first frame.  There is a weird',
     &                /,' bug causing multiple calls to <induce> to',
     &                /,' take a weirdly long amount of time.  You are',
     &                /,' better off separating the coordinate file',
     &                /,' into individual frame files.')
               goto 999
            end if
         end if
c         if (frame .ge. 10) then
c            print *, "Stopping at frame 10 for testing purposes"
c            goto 999
c         end if
c
c        rotate multipole components into global frame
c
         call rotpole
c
c        handle induced dipoles
c
         if (douind) then
            if (use_bounds .and. .not.use_rigid) then 
               call bounds
            end if
            if (use_list) call nblist
            call induce
c
c           if PB implicit solvent used, set uind to uinds
c
            if (solvtyp.eq.'PB') then
               do i = 1, npole
                  do j = 1, npole
                     uind(j,i) = uinds(j,i)
                  end do
               end do
            end if
            if (solvtyp.eq.'PB-HPMF') then
               do i = 1, npole
                  do j = 1, npole
                     uind(j,i) = uinds(j,i)
                  end do
               end do
            end if
c
c           Write to file to avoid recalculating
c
            lext = 3
            call numeral(frame, ext, lext)
            iind = freeunit()
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
            write(iind,200) n,title(1:ltitle)
  200       format(i0,2x,a)
            do i = 1, npole
               if (polarity(i) .ne. 0.0d0) then
                  k = ipole(i)
                  write (iind, 210) k,name(k),(debye*uind(j,i),j=1,3)
  210             format(i6,2x,a3,3f12.6)
               end if
            end do
            close(unit=iind)
            write(iout,220) indfile(1:trimtext(indfile))
  220       format(' Induced Dipole File', 10x, a)
         else 
c
c           read first line and return if already at eof
c
            quit = .true. 
            i = 0
            uindsize = 0
            do while (uindsize .eq. 0)
               read (iind,230,err=260,end=260) record
  230          format(a120)
               uindsize = trimtext(record)
            end do
            quit = .true. 
c
c           parse the title line to get the number of atoms
c
            next = 1
            call gettext(record,string,next)
            read(string,*,err=260,end=260) uindn
c
c           check for too many or too few atoms
c
            if (uindn .ne. n) then
               write(iout,240) n,uindn
  240          format(/,' Error, expected ',i0,' atoms in '
     &         'induced dipole file, but found ',i0)
               call fatal
            end if
c
c           initialize some variables
c
            do i = 1, uindn
               uindtag(i) = 0
               uindname(i) = ' '
            end do
c
c           read the coordinates
c        
            do i = 1, uindn
               uindsize = 0
               do while (uindsize .eq. 0)
                  read (iind,250,err=260,end=260) record
  250             format(a120)
                  uindsize = trimtext(record)
               end do
               read (record,*,err=260,end=260) uindtag(i)
               next = 1
               call getword(record,uindname(i),next)
               string = record(next:120)
               read(string,*,err=260,end=260) 
     &               uind(1,i),uind(2,i),uind(3,i)
               uind(1,i) = uind(1,i) / debye
               uind(2,i) = uind(2,i) / debye
               uind(3,i) = uind(3,i) / debye
            end do
            quit = .false.
  260       continue
            if (quit) then
               write(iout,270) i
  270          format(/,' Error in induced dipole file at Atom',i6)
               call fatal
            end if
         end if
         bondvector(1) = x(v2i) - x(v1i)
         bondvector(2) = y(v2i) - y(v1i)
         bondvector(3) = z(v2i) - z(v1i)
         bondlength = sqrt(bondvector(1)*bondvector(1)
     &                    +bondvector(2)*bondvector(2)
     &                    +bondvector(3)*bondvector(3))
c
c        Field = uind / polarizability
c
         call projectVector(uind(1:3,v1i),bondvector,v1field)
         v1field = v1field * kfac / polarity(v1i)
         call projectVector(uind(1:3,v2i),bondvector,v2field)
         v2field = v2field * kfac / polarity(v2i)
c
c        initialize some variables
c
         do i = 1, 9
            field(i) = 0.0d0
            exfield(i) = 0.0d0 
         end do
c
c        Find the field at a specific xyz coordinate due
c
         cx = (x(v1i)+x(v2i)) * 0.5d0
         cy = (y(v1i)+y(v2i)) * 0.5d0
         cz = (z(v1i)+z(v2i)) * 0.5d0
c
c        Find the field due to every atom
c
         do i = 1, n
            call dofield(i, cx, cy, cz, ifield, .false.)
            do j = 1, 9
               field(j) = field(j) + ifield(j)
            end do
         end do
c
c        Project along the bond vector
c
         call projectVector(field(1:3),bondvector,prfield(1))
         call projectVector(field(4:6),bondvector,prfield(2))
         call projectVector(field(7:9),bondvector,prfield(3))
         pfield = prfield(1) + prfield(2) + prfield(3)
         write(ifld,285,advance='no') frame,(v1field+v2field)*0.5,pfield
  285    format(' ',i8,2x,es16.6,2x,es16.6,2x)
c
c        Find the perm contributions of the excluded atoms
c
         do i = 1, nexatoms
            call dofield(exi(i), cx, cy, cz, ifield, .true.)
            do j = 1, 9
               exfield(j) = exfield(j) + ifield(j)
            end do
            call projectVector(ifield(1:3),bondvector,prexfield(1))
            call projectVector(ifield(4:6),bondvector,prexfield(2))
            call projectVector(ifield(7:9),bondvector,prexfield(3))
            write(ifld,286,advance='no') prexfield(1:3)
  286       format(es16.6,2x,es16.6,2x,es16.6,2x)
         end do
c
c        Project along the bond vector
c
         call projectVector(exfield(1:3),bondvector,prexfield(1))
         call projectVector(exfield(4:6),bondvector,prexfield(2))
         call projectVector(exfield(7:9),bondvector,prexfield(3))
         pexfield = prexfield(1) + prexfield(2) + prexfield(3)
c
c    Write to screen and file
c
c         write(iout,280) frame,(v1field+v2field)*0.5,pfield-pexfield,
c     &                   pfield
c  280    format(' ',i8,2x,es16.6,2x,es16.6,2x,es16.6)
         write(ifld,290,advance='no') pfield-pexfield
  290    format(es16.6,2x)
c
c        Find the total contributions of the excluded atoms
c
         do i = 1, 9
            exfield(i) = 0.0d0
         end do
         do i = 1, nexatoms
            call dofield(exi(i), cx, cy, cz, ifield, .false.)
            do j = 1, 9
               exfield(j) = exfield(j) + ifield(j)
            end do
         end do
c
c        Project along the bond vector
c
         call projectVector(exfield(1:3),bondvector,prexfield(1))
         call projectVector(exfield(4:6),bondvector,prexfield(2))
         call projectVector(exfield(7:9),bondvector,prexfield(3))
         pexfield = prexfield(1) + prexfield(2) + prexfield(3)
c
c    Write to screen and file
c
         write(ifld,295) pfield-pexfield
  295    format(es16.6)
c
c    attempt to read the next frame
c
         call readxyz(ixyz)
      end do
  999 continue
      close(unit=ixyz)
      close(unit=ifld)
      if (douind) close(unit=iind)
      end

      subroutine dofield(index,cx,cy,cz,field,uperm_only)
      use sizes
      use atoms
      use mpole
      use polar
      use units
      use chgpot
      use bound
      use boxes
      use cell
      integer index,i,j,k
      real*8 cx, cy, cz
      real*8 r,ri,r2,rr3,rr5,rr7,xr,yr,zr,pdotr,qdotr
      real*8 qq(3)
      real*8 cfac,coulombs,kbt,permitivity_constant
      real*8 field(9)
      logical uperm_only
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
      do i = 1, 3
         field(i) = 0.0d0
         field(i+3) = 0.0d0
         field(i+6) = 0.0d0
         qq(i) = 0.0d0
      end do
      xr = cx - x(index)
      yr = cy - y(index)
      zr = cz - z(index)
      r2 = xr*xr + yr*yr + zr*zr
      r = sqrt(r2)
      if (r2 .eq. 0.0d0) return
      if (use_bounds) call image(xr, yr, zr)
      r2 = xr*xr + yr*yr + zr*zr
      ri = sqrt(r2)
      if ((ri-r) .gt. 1.0d-6) then
            write(iout,10) index, r, ri
  10        format(' Error, atom ',i0,' should be <= ',f0.2,' Angstroms',
     &             /,' from point, but was found to be ',f0.2,
     &             /,' Angstroms away.')
            call fatal
      end if
      rr3 = 1.0d0 / (r*r2)
      rr5 = 3.0d0 / (r*r2*r2)
      rr7 = 15.0d0 / (r*r2*r2*r2)
c
c     Monopole field
c
      field(1) = rr3 * rpole(1,index) * xr
      field(2) = rr3 * rpole(1,index) * yr
      field(3) = rr3 * rpole(1,index) * zr
c
c     Dipole field
c
      if (.not. uperm_only) then
         pdotr = (rpole(2,index) + uind(1,index))*xr
     &         + (rpole(3,index) + uind(2,index))*yr
     &         + (rpole(4,index) + uind(3,index))*zr
         field(4) = rr5 * pdotr * xr - rr3
     &            * (rpole(2,index)+uind(1,index))
         field(5) = rr5 * pdotr * yr - rr3
     &            * (rpole(3,index)+uind(2,index))
         field(6) = rr5 * pdotr * zr - rr3
     &            * (rpole(4,index)+uind(3,index))
      else
         pdotr = (rpole(2,index))*xr
     &         + (rpole(3,index))*yr
     &         + (rpole(4,index))*zr
         field(4) = rr5 * pdotr * xr - rr3 * (rpole(2,index))
         field(5) = rr5 * pdotr * yr - rr3 * (rpole(3,index))
         field(6) = rr5 * pdotr * zr - rr3 * (rpole(4,index))
      end if
c
c     Quadrupole field
c
      qq(1) = rpole(5,index)*xr +
     &        rpole(6,index)*yr +
     &        rpole(7,index)*zr 
      qq(2) = rpole(6,index)*xr + 
     &        rpole(9,index)*yr +
     &        rpole(10,index)*zr 
      qq(3) = rpole(7,index)*xr + 
     &        rpole(10,index)*yr +
     &        rpole(13,index)*zr 
      qdotr = qq(1)*xr + qq(2)*yr + qq(3)*zr
      field(7) = rr7 * qdotr * xr - rr5 * qq(1) * 2.0d0
      field(8) = rr7 * qdotr * yr - rr5 * qq(2) * 2.0d0
      field(9) = rr7 * qdotr * zr - rr5 * qq(3) * 2.0d0
c
c     Total field = monopole + dipole + quadrupole
c
      do i = 1, 9
         field(i) = field(i) * cfac
      end do
      end 

      subroutine projectVector(field, vector, projection)
      real*8 field(3), vector(3)
      real*8 projection, vectorLength
      projection = 0.0d0
      vectorLength = sqrt(vector(1)*vector(1)
     &                  + vector(2)*vector(2)
     &                  + vector(3)*vector(3))
      projection = field(1)*vector(1)
     &           + field(2)*vector(2)
     &           + field(3)*vector(3)
      projection = projection / vectorLength
      end
