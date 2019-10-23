!
!-----------------------------------------------------------------------
!
!     auxillary routines that involve off-processor communication.
!     shared memory version, contained in module mod_xc.
!
!     author:  Alan J. Wallcraft,  NRL.
!
!-----------------------------------------------------------------------
!
      subroutine xcaget(aa, a, mnflg)
      implicit none
!
      real,    intent(out)   :: aa(itdm,jtdm)
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) convert an entire 2-D array from tiled to non-tiled layout.
!
!  3) mnflg selects which nodes must return the array
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  4) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aa              real           output    non-tiled target array
!    a               real           input     tiled source array
!    mnflg           integer        input     node return flag
!*
!**********
!
      integer j
#if defined(TIMER)
!
!     call xctmr0( 1)
#endif
!
!     use xclget for now.
!
      do j= 1,jtdm
        call xclget(aa(1,j),itdm, a, 1,j,1,0, mnflg)
      enddo
#if defined(TIMER)
!
!     call xctmr1( 1)
#endif
      return
      end subroutine xcaget

      subroutine xcaput(aa, a, mnflg)
      implicit none
!
      real,    intent(inout) :: aa(itdm,jtdm)
      real,    intent(out)   :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) convert an entire 2-D array from non-tiled to tiled layout.
!
!  3) mnflg selects which nodes must contain the non-tiled array
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  4) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aa              real           input     non-tiled source array
!    a               real           output    tiled target array
!    mnflg           integer        input     node source flag
!*
!**********
!
      integer j
#if defined(TIMER)
!
!     call xctmr0( 4)
#endif
!
!     use xclput for now.
!
      do j= 1,jtdm
        call xclput(aa(1,j),itdm, a, 1,j,1,0)
      enddo
#if defined(TIMER)
!
!     call xctmr1( 4)
#endif
      return
      end subroutine xcaput

      subroutine xcastr(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) broadcast array a to all tiles.
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    mnflg           integer        input     node originator flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0( 9)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1( 9)
#endif
      return
      end subroutine xcastr

      subroutine xceget(aelem, a, ia,ja, mnflg)
      implicit none
!
      real,    intent(out)   ::         aelem
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    ::         ia,ja
      integer, intent(in), optional ::  mnflg
!
!**********
!*
!  1) find the value of a(ia,ja) on the non-tiled 2-D grid.
!
!  2) mnflg selects which nodes must return the line
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aelem           real           output    required element
!    a               real           input     source array
!    ia              integer        input     1st index into a
!    ja              integer        input     2nd index into a
!    mnflg           integer        input     node return flag
!*
!**********
#if defined(TIMER)
!
      call xctmr0( 2)
#endif
!
!     single node version - trivial indexing.
!
      aelem = a(ia,ja)
#if defined(TIMER)
!
      call xctmr1( 2)
#endif
      return
      end subroutine xceget

      subroutine xceput(aelem, a, ia,ja, mnflg)
      implicit none
!
      integer, intent(in), optional ::  mnflg
      integer, intent(in)    ::         ia,ja
      real,    intent(in)    ::         aelem
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) fill a single element in the non-tiled 2-D grid.
!
!  2) mnflg selects which nodes hold the element on entry
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aelem           real           input     element value
!    a               real           in/out    target array
!    ia              integer        input     1st index into a
!    ja              integer        input     2nd index into a
!    mnflg           integer        input     node source flag
!*
!**********
#if defined(TIMER)
!
      call xctmr0( 4)
#endif
!
!     single node version - trivial indexing.
!
      a(ia,ja) = aelem
#if defined(TIMER)
!
      call xctmr1( 4)
#endif
      return
      end subroutine xceput

      subroutine xchalt(cerror)
      implicit none
!
      character*(*), intent(in) :: cerror
!
!**********
!*
!  1) stop all processes.
!
!  2) only one processes need call this routine, i.e. it is for
!      emergency stops.  use 'xcstop' for ordinary stops called
!      by all processes.
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    cerror          char*(*)       input     error message
!*
!**********
!
!     shared memory version, just stop.
!
      if     (cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
        write(lp,*) '**************************************************'
        call flush(lp)
      endif
      stop '(xchalt)'
      end subroutine xchalt

      subroutine xciget(alist,nl, a, ia,ja, mnflg)
      implicit none
!
      integer, intent(in)    ::         nl
      real,    intent(out)   ::         alist(nl)
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    ::         ia(nl),ja(nl)
      integer, intent(in), optional ::  mnflg
!
!**********
!*
!  1) find the value of a(ia(:),ja(:)) on the non-tiled 2-D grid.
!
!  2) mnflg selects which nodes must return the list of elements
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    alist           real           output    required elements
!    nl              integer        input     dimension of alist
!    a               real           input     source array
!    ia              integer        input     1st indexes into a
!    ja              integer        input     2nd indexes into a
!    mnflg           integer        input     node return flag
!*
!**********
!
      integer i
!
#if defined(TIMER)
!
      call xctmr0( 2)
#endif
!
!     single node version - trivial indexing.
!
      do i= 1,nl
        alist(i) = a(ia(i),ja(i))
      enddo !i
#if defined(TIMER)
!
      call xctmr1( 2)
#endif
      return
      end subroutine xciget

      subroutine xciget_sm(alist,nl, a, ia,ja, mnflg)
      implicit none
!
      integer, intent(in)    ::         nl
      real,    intent(out)   ::         alist(nl)
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    ::         ia(nl),ja(nl)
      integer, intent(in), optional ::  mnflg
!
!**********
!*
!  1) find the value of a(ia(:),ja(:)) on the non-tiled 2-D grid.
!     identical to xciget and never invoked (for compatibility only)
!
!  2) mnflg selects which nodes must return the list of elements
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    alist           real           output    required elements
!    nl              integer        input     dimension of alist
!    a               real           input     source array
!    ia              integer        input     1st indexes into a
!    ja              integer        input     2nd indexes into a
!    mnflg           integer        input     node return flag
!*
!**********
!
      integer i
!
#if defined(TIMER)
!
      call xctmr0( 2)
#endif
!
!     single node version - trivial indexing.
!
      do i= 1,nl
        alist(i) = a(ia(i),ja(i))
      enddo !i
#if defined(TIMER)
!
      call xctmr1( 2)
#endif
      return
      end subroutine xciget_sm

      subroutine xciput(alist,nl, a, ia,ja, mnflg)
      implicit none
!
      integer, intent(in)    ::         nl
      integer, intent(in), optional ::  mnflg
      integer, intent(in)    ::         ia(nl),ja(nl)
      real,    intent(inout) ::         alist(nl)
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) fill a list of elements in the non-tiled 2-D grid.
!
!  2) mnflg selects which nodes hold the elements on entry
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    alist           real           in/out    element values
!    nl              integer        input     dimension of alist
!    a               real           in/out    target array
!    ia              integer        input     1st indexes into a
!    ja              integer        input     2nd indexes into a
!    mnflg           integer        input     node source flag
!*
!**********
!
      integer i
#if defined(TIMER)
!
      call xctmr0( 4)
#endif
!
!     single node version - trivial indexing.
!
      do i= 1,nl
        a(ia(i),ja(i)) = alist(i)
      enddo !1
#if defined(TIMER)
!
      call xctmr1( 4)
#endif
      return
      end subroutine xciput

      subroutine xclget(aline,nl, a, i1,j1,iinc,jinc, mnflg)
      implicit none
!
      integer, intent(in)    ::  nl,i1,j1,iinc,jinc,mnflg
      real,    intent(out)   ::  aline(nl)
      real,    intent(in)    ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) extract a line of elements from the non-tiled 2-D grid.
!
!  2) aline(i) = a(i1+iinc*(i-1),j1+jinc*(i-1)), for i=1...nl.
!     iinc and jinc can each be -1, 0, or +1.
!
!     if jinc=0, j1 can be between jtdm+1 and jtdm+nbdy to return
!     values from the top halo.  This is for debugging the arctic
!     patch halo exchange only.
!
!  3) mnflg selects which nodes must return the line
!        =-n; node number n (mnproc=n), nl,i1,j1 only on node n
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!     normally all integer arguments must be identical on all nodes,
!     but a negative mnflg indicates that only the target node is
!     providing the nl,i1,j1 values.  These are broadcast to all other
!     nodes, and returned in nl,i1,j1 by all of them.
!
!     mnflg is ignored here (only have a single node).
!
!  4) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aline           real           output    required line of elements
!    nl              integer        in/out    dimension of aline
!    a               real           input     source array
!    i1              integer        in/out    1st index into a
!    j1              integer        in/out    2nd index into a
!    iinc            integer        input     1st index increment
!    jinc            integer        input     2nd index increment
!    mnflg           integer        input     node return flag
!
!    nl,i1,j1 are input only unless mnflg is negative.
!*
!**********
!
      integer i
#if defined(TIMER)
!
      call xctmr0( 3)
#endif
!
!     single node version - trivial indexing and no error checking.
!
      if     (jinc.eq.0) then
        do i= 1,nl
          aline(i) = a(i1+iinc*(i-1),j1)
        enddo
      elseif (iinc.eq.0) then
        do i= 1,nl
          aline(i) = a(i1,j1+jinc*(i-1))
        enddo
      else
        do i= 1,nl
          aline(i) = a(i1+iinc*(i-1),j1+jinc*(i-1))
        enddo
      endif
#if defined(TIMER)
!
      call xctmr1( 3)
#endif
      return
      end subroutine xclget

      subroutine xclput(aline,nl, a, i1,j1,iinc,jinc, mnflg)
      implicit none
!
      integer, intent(in), optional ::  mnflg
      integer, intent(in)    ::         nl,i1,j1,iinc,jinc
      real,    intent(inout) ::         aline(nl)
      real,    intent(inout) ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) fill a line of elements in the non-tiled 2-D grid.
!
!  2) aline(i) = a(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
!     one of iinc and jinc must be 0, and the other must be 1.
!
!  3) mnflg selects which nodes hold the line on entry
!        = 0; all nodes (default)
!        = n; node number n (mnproc=n)
!     mnflg is ignored here (only have a single node).
!
!  4) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    aline           real           in/out    line of element values
!    nl              integer        input     dimension of aline
!    a               real           in/out    target array
!    i1              integer        input     1st index into a
!    j1              integer        input     2nd index into a
!    iinc            integer        input     1st index increment
!    jinc            integer        input     2nd index increment
!    mnflg           integer        input     node source flag
!*
!**********
!
      integer i
#if defined(TIMER)
!
      call xctmr0( 4)
#endif
!
!     single node version - trivial indexing.
!
      if     (jinc.eq.0) then
        do i= 1,nl
          a(i1+i-1,j1) = aline(i)
        enddo
      elseif (iinc.eq.0) then
        do i= 1,nl
          a(i1,j1+i-1) = aline(i)
        enddo
      endif
#if defined(TIMER)
!
      call xctmr1( 4)
#endif
      return
      end subroutine xclput

      subroutine xcmaxr_0(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace scalar a with its element-wise maximum over all tiles.
!
!  2) mnflg selects which nodes must return the minimum
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target variable
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcmaxr_0

      subroutine xcmaxr_1(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace array a with its element-wise maximum over all tiles.
!
!  2) mnflg selects which nodes must return the minimum
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcmaxr_1

      subroutine xcmaxr_0o(a)
      implicit none
!
      real, intent(inout) :: a
!
!**********
!*
!  1) replace scalar a with its element-wise maximum over all tiles.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target variable
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcmaxr_0o

      subroutine xcmaxr_1o(a)
      implicit none
!
      real, intent(inout) :: a(:)
!
!**********
!*
!  1) replace array a with its element-wise maximum over all tiles.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcmaxr_1o

      subroutine xcminr_0(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace scalar a with its element-wise minimum over all tiles.
!
!  2) mnflg selects which nodes must return the minimum
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target variable
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcminr_0

      subroutine xcminr_1(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace array a with its element-wise minimum over all tiles.
!
!  2) mnflg selects which nodes must return the minimum
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcminr_1

      subroutine xcminr_0o(a)
      implicit none
!
      real, intent(inout) :: a
!
!**********
!*
!  1) replace scalar a with its element-wise minimum over all tiles.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target variable
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcminr_0o

      subroutine xcminr_1o(a)
      implicit none
!
      real, intent(inout) :: a(:)
!
!**********
!*
!  1) replace array a with its element-wise minimum over all tiles.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcminr_1o

      subroutine xcspmd
      implicit none
!
!**********
!*
!  1) initialize data structures that identify the tiles.
!
!  2) data structures:
!      ipr     - 1st 2-D node dimension
!      jpr     - 2nd 2-D node dimension
!      ijpr    -     1-D node dimension (ipr*jpr)
!      mproc   - 1st 2-D node index
!      nproc   - 2nd 2-D node index
!      mnproc  -     1-D node index
!      i0      -     1st dimension tile offset
!      ii      -     1st dimension tile extent
!      j0      -     2nd dimension tile offset
!      jj      -     2nd dimension tile extent
!      nreg    -     region type
!      vland   -     fill value for land (standard value 0.0)
!
!  3) ipr,jpr,ijpr are global (tile independent) values.
!     all other values depend on the processor number, 
!     but in this case there is only one processor.
!*
!**********
!
!     shared memory version, mproc=nproc=1.
!
      if     (iqr.ne.1 .or. jqr.ne.1 .or. ijqr.ne.1) then
        call xcstop('Error in xcspmd: must have iqr=jqr=ijqr=1')
               stop '(xcspmd)'
      endif
!
      lp = 6
!
      ipr    = 1
      jpr    = 1
      ijpr   = 1
      mnproc = 1
      mproc  = 1
      nproc  = 1
!
#if defined(RELO)
! --- region's horizontal dimensions are from blkdat.input.
!
      itdm = -1
      jtdm = -1
#endif
!
      i0  = 0
      ii  = itdm
      j0  = 0
      jj  = jtdm
!
      nreg   = -1  ! unknown region type
!
      vland  = 0.0
      vland4 = 0.0
!
!     initialize timers.
!
      call xctmri
#if defined(TIMER)
      call xctmrn( 1,'xcaget')
      call xctmrn( 2,'xceget')
      call xctmrn( 3,'xclget')
      call xctmrn( 4,'xcXput')
      call xctmrn( 5,'xcsum ')
      call xctmrn( 6,'xcrang')
      call xctmrn( 9,'xcastr')
      call xctmrn(10,'xcmaxr')
      call xctmrn(12,'xctilr')
#if defined(ARCTIC)
      call xctmrn(13,'xctila')
#endif
#endif
      return
      end subroutine xcspmd

      subroutine xcstop(cerror)
      implicit none
!
      character*(*), intent(in) :: cerror
!
!**********
!*
!  1) stop all processes.
!
!  2) all processes must call this routine.
!     use 'xchalt' for emergency stops.
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    cerror          char*(*)       input     error message
!*
!**********
!
!     print active timers.
!
      call xctmrp
!
!     shared memory version, just stop.
!
      if     (cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
        write(lp,*) '**************************************************'
        call flush(lp)
      endif
#if defined(STOP2003)
! --- Fortran 2003 STOP outputs many lines, replace with exit
      call exit(0)
#else
      stop '(xcstop)' 
#endif
      end subroutine xcstop

      subroutine xcsum(sum, a,mask)
      implicit none
!
      real*8,  intent(out)   :: sum
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) sum a 2-d array, where mask==1
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    sum             real*8         output    sum of a
!    a               real           input     source array
!    mask            integer        input     mask array
!
!  3) sum is bit for bit reproducable for the same halo size, nbdy.
!*
!**********
!
      real*8     zero8
      parameter (zero8=0.0)
!
      real*8  sum8,sum8p,sum8j(jdm)
      integer i,i1,j
#if defined(TIMER)
!
      call xctmr0( 5)
#endif
!
!     row sums in 2*nbdy+1 wide strips.
!
!$OMP PARALLEL DO PRIVATE(j,i1,i,sum8,sum8p) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jdm
        sum8 = zero8
        do i1=1,idm,2*nbdy+1
          sum8p = zero8
          do i= i1,min(i1+2*nbdy,idm)
            if     (mask(i,j).eq.1) then
              sum8p = sum8p + a(i,j)
            endif
          enddo
          sum8 = sum8 + sum8p
        enddo
        sum8j(j) = sum8  ! use of sum8 minimizes false sharing of sum8j
      enddo
!$OMP END PARALLEL DO
!
!     serial sum of rwo-sum loop.
!
      sum8 = sum8j(1)
      do j=2,jdm
        sum8 = sum8 + sum8j(j)
      enddo
      sum = sum8
#if defined(TIMER)
!
      call xctmr1( 5)
#endif
      return
      end subroutine xcsum

      subroutine xcsumj(sumj, a,mask)
      implicit none
!
      real*8,  intent(out)   :: sumj(jtdm)
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
!**********
!*
!  1) rwo-sum of a 2-d array, where mask==1, on first processor only
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    sumj            real*8         output    row-sum of a
!    a               real           input     source array
!    mask            integer        input     mask array
!
!  3) sum is bit for bit reproducable for the same halo size, nbdy.
!*
!**********
!
      real*8     zero8
      parameter (zero8=0.0)
!
      real*8  sum8,sum8p
      integer i,i1,j
#if defined(TIMER)
!
      call xctmr0( 5)
#endif
!
!     row sums in 2*nbdy+1 wide strips.
!
!$OMP PARALLEL DO PRIVATE(j,i1,i,sum8,sum8p) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jdm
        sum8 = zero8
        do i1=1,idm,2*nbdy+1
          sum8p = zero8
          do i= i1,min(i1+2*nbdy,idm)
            if     (mask(i,j).eq.1) then
              sum8p = sum8p + a(i,j)
            endif
          enddo
          sum8 = sum8 + sum8p
        enddo
        sumj(j) = sum8  ! use of sum8 minimizes false sharing of sumj
      enddo
!$OMP END PARALLEL DO
#if defined(TIMER)
!
      call xctmr1( 5)
#endif
      return
      end subroutine xcsumj

      subroutine xcsumr_0(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace scalar a with its element-wise sum over all tiles.
!
!  2) mnflg selects which nodes must return the sum 
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target variable
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcsumr_0

      subroutine xcsumr_1(a, mnflg)
      implicit none
!
      real,    intent(inout) :: a(:)
      integer, intent(in)    :: mnflg
!
!**********
!*
!  1) replace array a with its element-wise sum over all tiles.
!
!  2) mnflg selects which nodes must return the sum 
!        = 0; all nodes
!        = n; node number n (mnproc=n)
!
!  3) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    mnflg           integer        input     node return flag
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(10)
#endif
!
!     single node version - do nothing.
#if defined(TIMER)
!
      call xctmr1(10)
#endif
      return
      end subroutine xcsumr_1

      subroutine xcsync(lflush)
      implicit none
!
      logical, intent(in) :: lflush
!
!**********
!*
!  1) barrier, no processor exits until all arrive (and flush stdout).
!
!  2) some MPI implementations only flush stdout as a collective
!     operation, and hence the lflush=.true. option to flush stdout.
!
!  3) Only one processor, so the barrier is a no-op in this case.
!*
!**********
!
      if     (lflush) then
        call flush(lp)
      endif
      return
      end subroutine xcsync

#if defined(ARCTIC)
      subroutine xctila(a,l1,ld,itype)
      implicit none
!
      integer, intent(in)    :: l1,ld,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
!
!**********
!*
!  1) update the top row of a real array.
!     only needed when importing a tripole grid array.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    l1              integer        input     3rd dim. start index
!    ld              integer        input     3rd dimension of a
!    itype           integer        input     grid and field type
!
!  3) itype selects both the grid and field type
!        itype= 1; p-grid, scalar field
!        itype=11; p-grid, vector field
!
!  4) this version for a global tripole (arctic bipolar patch) grid
!*
!**********
!
      integer i,io,k
#if defined(TIMER)
!
      call xctmr0(13)
#endif
      do k= l1,ld
        if     (itype.eq.1) then
!
!         scalar on p-grid
!
          do i= 1,ii
            io = ii-mod(i-1,ii)
            a(i,jj,k) = a(io,jj-1,k)
          enddo
        else
!
!         vector p-grid field, swap sign
!
          do i= 1,ii
            io = ii-mod(i-1,ii)
            a(i,jj,k) = -a(io,jj-1,k)
          enddo
        endif
      enddo
#if defined(TIMER)
!
      call xctmr1(13)
#endif
      return
      end subroutine xctila
      subroutine xctilr(a,l1,ld,mh,nh,itype)
      implicit none
!
      integer, intent(in)    :: l1,ld,mh,nh,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
!
!**********
!*
!  1) update the tile overlap halo of a real array.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    l1              integer        input     3rd dim. start index
!    ld              integer        input     3rd dimension of a
!    mh              integer        input     1st (EW) update halo size
!    nh              integer        input     2nd (NS) update halo size
!    itype           integer        input     grid and field type
!
!  3) itype selects both the grid and field type
!        itype= 1; p-grid, scalar field
!        itype= 2; q-grid, scalar field
!        itype= 3; u-grid, scalar field
!        itype= 4; v-grid, scalar field
!        itype=11; p-grid, vector field
!        itype=12; q-grid, vector field
!        itype=13; u-grid, vector field
!        itype=14; v-grid, vector field
!
!  4) this version for a global grid that includes the arctic ocean
!*
!**********
!
      integer i,io,j,k,mhl,nhl
#if defined(TIMER)
!
      call xctmr0(12)
#endif
!
      mhl = max(0,min(mh,nbdy))
      nhl = max(0,min(nh,nbdy))
!
      if     (nhl.gt.0) then
        do k= l1,ld
!
!         southern boundary is closed.
!
          do j= 1,nhl
            do i= 1,ii
              a(i,1-j,k) = vland
            enddo
          enddo
!
          if     (itype.lt.10) then
!
!           scalar field
!
            if     (itype.eq. 1) then
!
!             p-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = ii-mod(i-1,ii)
                  a(i,jj+j,k) = a(io,jj-1-j,k)
                enddo
              enddo
            elseif (itype.eq. 2) then
!
!             q-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = mod(ii-(i-1),ii)+1
                  a(i,jj+j,k) = a(io,jj-j,k)
                enddo
              enddo
            elseif (itype.eq. 3) then
!
!             u-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = mod(ii-(i-1),ii)+1
                  a(i,jj+j,k) = a(io,jj-1-j,k)
                enddo
              enddo
            else
!
!             v-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = ii-mod(i-1,ii)
                  a(i,jj+j,k) = a(io,jj-j,k)
                enddo
              enddo
            endif
          else
!
!           vector field, swap sign
!
            if     (itype.eq.11) then
!
!             p-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = ii-mod(i-1,ii)
                  a(i,jj+j,k) = -a(io,jj-1-j,k)
                enddo
              enddo
            elseif (itype.eq.12) then
!
!             q-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = mod(ii-(i-1),ii)+1
                  a(i,jj+j,k) = -a(io,jj-j,k)
                enddo
              enddo
            elseif (itype.eq.13) then
!
!             u-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = mod(ii-(i-1),ii)+1
                  a(i,jj+j,k) = -a(io,jj-1-j,k)
                enddo
              enddo
            else
!
!             v-grid
!
              do j= 1,nhl
                do i= 1,ii
                  io = ii-mod(i-1,ii)
                  a(i,jj+j,k) = -a(io,jj-j,k)
                enddo
              enddo
            endif
          endif
        enddo
      endif
!
      if     (mhl.gt.0) then
        do k= 1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              a( 1-i,j,k) = a(ii+1-i,j,k)
              a(ii+i,j,k) = a(     i,j,k)
            enddo
          enddo
        enddo
      endif
#if defined(TIMER)
!
      call xctmr1(12)
#endif
      return
      end subroutine xctilr
#else /* !ARCTIC */
      subroutine xctilr(a,l1,ld,mh,nh,itype)
      implicit none
!
      integer, intent(in)    :: l1,ld,mh,nh,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
!
!**********
!*
!  1) update the tile overlap halo of a real   array.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    a               real           in/out    target array
!    l1              integer        input     3rd dim. start index
!    ld              integer        input     3rd dimension of a
!    mh              integer        input     1st (EW) update halo size
!    nh              integer        input     2nd (NS) update halo size
!    itype           integer        input     grid and field type
!
!  3) itype selects both the grid and field type
!        itype= 1; p-grid, scalar field
!        itype= 2; q-grid, scalar field
!        itype= 3; u-grid, scalar field
!        itype= 4; v-grid, scalar field
!        itype=11; p-grid, vector field
!        itype=12; q-grid, vector field
!        itype=13; u-grid, vector field
!        itype=14; v-grid, vector field
!     it is ignored here because all types are the same unless
!      the grid includes the arctic ocean
!*
!**********
!
      integer i,j,k,mhl,nhl
#if defined(TIMER)
!
      call xctmr0(12)
#endif
!
      mhl = max(0,min(mh,nbdy))
      nhl = max(0,min(nh,nbdy))
!
      if     (nhl.gt.0) then
        if     (nreg.le.2) then  ! closed in latitude
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = vland
                a(i,jj+j,k) = vland
              enddo
            enddo
          enddo
        else  ! periodic (f-plane) in latitude
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = a(i,jj+1-j,k)
                a(i,jj+j,k) = a(i,     j,k)
              enddo
            enddo
          enddo
        endif
      endif
!
      if     (mhl.gt.0) then
        if     (nreg.eq.0 .or. nreg.eq.4) then  ! closed in longitude
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = vland
                a(ii+i,j,k) = vland
              enddo
            enddo
          enddo
        else  ! periodic in longitude
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = a(ii+1-i,j,k)
                a(ii+i,j,k) = a(     i,j,k)
              enddo
            enddo
          enddo
        endif
      endif
#if defined(TIMER)
!
      call xctmr1(12)
#endif
      return
      end subroutine xctilr
#endif /* ARCTIC:else */

      subroutine xctmri
      implicit none
!
!**********
!*
!  1) initialize timers.
!
!  2) timers  1:32 are for message passing routines,
!     timers 33:80 are for general hycom routines,
!     timers 81:96 are for user selected routines.
!     timer     97 is the total time.
!
!  3) call xctmri    to initialize timers (called in xcspmd),
!     call xctmr0(n) to start timer n,
!     call xctmr1(n) to stop  timer n and add event to timer sum,
!     call xctnrn(n,cname) to register a name for timer n,
!     call xctmrp to printout timer statistics (called by xcstop).
!*
!**********
!
      integer i
!
      real*8     zero8
      parameter (zero8=0.0)
!
      do 110 i= 1,97
        cc(i) = '      '
        nc(i) = 0
        tc(i) = zero8
  110 continue
!
      call xctmrn(97,'total ')
      call xctmr0(97)
      return
      end subroutine xctmri

      subroutine xctmr0(n)
      implicit none
!
      integer, intent(in) :: n
!
!**********
!*
!  1) start timer n.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    n               integer        input     timer number
!
!  3) time every 50-th event above 5,000.
!*
!**********
!
      real*8 wtime
!
#if defined(DEBUG_TIMER)
      if     (n.gt.24 .and. cc(n).ne.'      ') then
        write(lp,*) 'call ',cc(n)
        call flush(lp)
      endif
#endif
      if     (timer_on) then
        if     (mod(nc(n),50).eq.0 .or. nc(n).le.5000) then
          t0(n) = wtime()
        endif
      endif !timer_on
      return
      end subroutine xctmr0

      subroutine xctmr1(n)
      implicit none
!
      integer, intent(in) :: n
!
!**********
!*
!  1) add time since call to xctim0 to timer n.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    n               integer        input     timer number
!                                                                        
!  3) time every 50-th event above 5,000.                 
!*
!**********
!
      real*8  wtime
!
      if     (timer_on) then
        if     (nc(n).gt.5000) then
          if     (mod(nc(n),50).eq.0) then
            tc(n) = tc(n) + 50.0*(wtime() - t0(n))
          endif                                   
        else                                      
          tc(n) = tc(n) + (wtime() - t0(n))
        endif                              
        nc(n) = nc(n) + 1                  
      endif !timer_on    
#if defined(DEBUG_TIMER)
      if     (n.gt.24 .and. cc(n).ne.'      ') then
        write(lp,*) 'exit ',cc(n)
        call flush(lp)
      endif
#endif
      return
      end subroutine xctmr1

      subroutine xctmrn(n,cname)
      implicit none
!
      character*6, intent(in) :: cname
      integer,     intent(in) :: n
!
!**********
!*
!  1) register name of timer n.
!
!  2) parameters:
!       name            type         usage            description
!    ----------      ----------     -------  ----------------------------
!    n               integer        input     timer number
!    cname           char*(8)       input     timer name
!*
!**********
!
      cc(n) = cname
      return
      end subroutine xctmrn

      subroutine xctmrp
      implicit none
!
!**********
!*
!  1) print all active timers.
!
!  2) on exit all timers are reset to zero.
!*
!**********
!
      real*8  tci
      integer i
!
      real*8     zero8
      parameter (zero8=0.0)
!
!     get total time.
!
      call xctmr1(97)
!
!     print timers.
!
      write(lp,6000)
      do i= 1,97
        if     (nc(i).ne.0) then
          if     (nc(i).le.5000) then
            tci = tc(i)
          else !correct for timing every 50th event
            tci = (tc(i)*nc(i))/(nc(i)-mod(nc(i),50))
          endif
          if     (cc(i).ne.'      ') then
            write(lp,6100) cc(i),nc(i),tci,tci/nc(i)
          else
            write(lp,6150)    i, nc(i),tci,tci/nc(i)
          endif
        endif
      enddo
      write(lp,6200)
      call flush(lp)
!
!     reset timers to zero.
!
      do i= 1,97
        nc(i) = 0
        tc(i) = zero8
      enddo
!
!     start a new total time measurement.
!
      call xctmr0(97)
      return
!
 6000 format(/ / &
          4x,' timer statistics ' / &
          4x,'------------------' /)
 6100 format(5x,a6, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
 6150 format(5x,'   #',i2, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
 6200 format(/ /)
      end subroutine xctmrp
!
!  Revision history:
!
!> Nov. 2011 - time every 50-th event above 5,000 (was 1,000).
!> Mar. 2012 - added optional mnflg to xclput
!> Apr. 2012 - added optional mnflg to xceget and xceput
!> Apr. 2012 - added xciget and xciput
!> Oct. 2019 - added xcsumr
