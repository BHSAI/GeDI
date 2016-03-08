      subroutine invrs(mat,N,nb,nproc,myid)
!-------------------------
!     Inverts N x N matrix 
!     mat : on input, matrix to be inverted; 1D array, row-major
!           on output, the inverse; 1D, row-major
!     N   : dimension
!     nb  : block size
!     nproc : number of processes
!     myid :  process id (0,1,...,nproc-1)
!     Adapted from
!       wwwuser.gwdg.de/~ohaan/ScaLAPACK_examples/use_PDGESV.f
!-------------------------
      implicit none 
      integer   N,nb
      real*8    mat(*)
      integer   nproc,myid

!  parameter constants
      integer    max_al, max_bl, max_vl
      parameter(max_al=10000000, max_bl=10000000,  max_vl=10000000)
!  variables for definition of matrix A and matrix X
!  local variables
      real*8    al(max_al), bl(max_bl), rl(max_bl)
      real*8    ai(max_bl),ainv(max_bl)
      integer   ctxt, ctxt_sys, ctxt_all, 
     :          nprow, npcol, myrow, mycol,
     :          m_al, n_al, m_bl, n_bl,
     :          llda, lldb, lldc, desc_A(9), desc_B(9),
     :          ncheck, k,l,
     :          irhs, ipr, ipc, il, jl, i, j, info, vl(max_vl),  
     :          NUMROC, INDXL2G, INDXG2L, INDXG2P
      call BLACS_GET( 0, 0, ctxt_sys )
      ctxt_all = ctxt_sys
      call BLACS_GRIDINIT( ctxt_all, 'C', nproc, 1)

! Parameters for processor grid
      nprow=int(sqrt(real(nproc)))
      npcol=int(real(nproc)/nprow)
! Set up a process grid of size nprow*npcol
      if (nprow*npcol.gt.nproc) then
        write(6,*) 'nproc = ',nproc,' less then nprow*npcol = '
        call BLACS_EXIT(1); stop
      end if
      ctxt = ctxt_sys; call BLACS_GRIDINIT( ctxt, 'C', nprow, npcol)

! Processes not belonging to the grid jump to the end of program
      if (ctxt.lt.0) go to 1000

! Get the process coordinates in the grid
      call BLACS_GRIDINFO( ctxt, nprow, npcol, myrow, mycol )
! number of rows and columns of local parts al, bl for A, B
      m_al = NUMROC( n, nb, myrow, 0, nprow )
      n_al = NUMROC( n, nb, mycol, 0, npcol )
      m_bl = NUMROC( N, nb, myrow, 0, nprow )
      n_bl = NUMROC( n, nb, mycol, 0, npcol )
!     write(*,*)myid,nprow,npcol
!     write(*,*)myid,myrow,mycol,m_al,n_al,m_bl,n_bl
!     call BLACS_EXIT(1)
!     stop
! Test for sufficient memory
      if (m_al*n_al.gt.max_al.or.m_bl*n_bl.gt.max_bl) then
        write(6,*)'not enough memory:  al ',m_al*n_al,max_al,
     :                                'bl ',m_bl*n_bl,max_bl
        call BLACS_EXIT(1)
        stop
      end if 

! initializing descriptors for the distributed matrices A, B:
      llda = max(1,m_al); lldb = max(1,m_bl)
      call DESCINIT( desc_A, n, n, nb, nb,0,0, ctxt, llda, info )
      call DESCINIT( desc_B, n, n, nb, nb,0,0, ctxt, lldb, info )

! initialize in parallel the local parts of A in al and in al_s
      do jl = 1 , n_al
         j = INDXL2G(jl,nb,mycol,0,npcol)
         do il = 1 , m_al
            i = INDXL2G(il,nb,myrow,0,nprow)
            al(il+m_al*(jl-1))= mat((i-1)*N+j)
!           write(*,*)j,i,mat((i-1)*N+j)
         end do
      end do
! initialize in parallel the local parts of B=I in bl
      do jl = 1 , n_bl
         j = INDXL2G(jl,nb,mycol,0,npcol)
         do il = 1 , m_bl
            i = INDXL2G(il,nb,myrow,0,nprow)
            if(i.eq.j) then
              bl(il+m_bl*(jl-1))=1
            else
              bl(il+m_bl*(jl-1))=0
            endif
         end do
      end do

! solve the systems of equations
      call PDGESV(n,n,al,1,1,desc_A,vl,bl,1,1,desc_B,info)

!     call BLACS_EXIT(1); stop

      if(myrow.eq.0 .and. mycol.eq.0) then
        do i=0,nprow-1
          do j=0,npcol-1
            m_bl = NUMROC(n,nb,i,0,nprow)
            n_bl = NUMROC(n,nb,j,0,npcol)
            if(i.ne.0 .or. j.ne.0) then
              call dgerv2d(ctxt,m_bl,n_bl,ai,m_bl,i,j)
!             write(*,*)'recvd',myid
            endif
            do jl=1,n_bl
              l = INDXL2G(jl,nb,j,0,npcol)
              do il=1,m_bl
                k = INDXL2G(il,nb,i,0,nprow)
                if(i.eq.0.and.j.eq.0) then
                  ainv(l+(k-1)*n)=bl(il+m_bl*(jl-1))
                else
                  ainv(l+(k-1)*n)=ai(il+m_bl*(jl-1))
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do jl=1,n_bl
          l = INDXL2G(jl,nb,mycol,0,npcol)
          do il=1,m_bl
            k = INDXL2G(il,nb,myrow,0,nprow)
          enddo
        enddo
        call dgesd2d(ctxt,m_bl,n_bl,bl,m_bl,0,0)
!       write(*,*)'sent',myid
      endif

      do j=1,n
        do i=1,n
          mat((j-1)*N+i)=ainv((j-1)*n+i)
        enddo
      enddo

 1000 continue
      end

