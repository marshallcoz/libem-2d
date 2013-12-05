!     Variables   
      module gloVars 
      save
      ! verbose = 0   ! no output 
      !         = 1   ! main calls, properties and images
      !         = 2   ! 1 + counters in loops and subfunctions
      !         = 3   ! 2 + matrix values
      integer    :: verbose 
      logical    :: makeVideo
      logical    :: workBoundary
      logical    :: plotFKS
      integer    :: multSubdiv ! Lo ideal es 4 o un multiplo
      logical    :: symmetricProblem
      real       :: periodicdamper
      logical, parameter    :: getInquirePointsSol = .true.!dont change
      integer, parameter    :: imecMax = 2
      real, dimension(2), parameter :: fixedNormal = (/0.0,1.0/)
      complex*16, parameter :: UI = cmplx(0.0,1.0,8), &
                               UR = cmplx(1.0,0.0,8)
      real*8, parameter :: PI = real(4.0*ATAN(1.0),8)
      integer, parameter :: Hplot = 700 , Wplot = 1200 
      end module gloVars
      
      module Gquadrature
      real :: WLmulti ! una o dos longitudes de onda
      integer, parameter :: Gquad_n = 8 ! número par
      
!     real, parameter, dimension(1) :: Gqu_t_8 = (/1.0/)
!     real, parameter, dimension(1) :: Gqu_A_8 =(/1.0/)
      real, parameter, dimension(8) :: Gqu_t_8 = & 
      (/ -0.96028986, &
         -0.79666648, &
         -0.52553242, &
         -0.18343464, &
         +0.18343464, &
         +0.52553242, &
         +0.79666648, &
         +0.96028986 /)
      
      real, parameter, dimension(8) :: Gqu_A_8 = &
      (/ 0.10122854, &
         0.22238104, &
         0.31370664, &
         0.36268378, &
         0.36268378, &
         0.31370664, &
         0.22238104, &
         0.10122854 /)
         
      end module Gquadrature 
      
      module soilVars
!     use gloVars, only: dp
      save
      integer ::  N !number of layers. HALF-SPACE at N+1
      real, dimension(:),  allocatable :: Z,AMU,BETA,ALFA,LAMBDA,RHO               
      end module soilVars
      
      module sourceVars
      integer :: efsource
      real    :: xfsource,zfsource,nxfsource,nzfsource
      logical :: intfsource
      end module sourceVars
      
      module waveNumVars
      !frequency loop vars:
      integer,save      :: NFREC
      real   ,save      :: FREC,DFREC,OME,OMEI
      !Discrete Wave-number:
      real   ,save      :: DK    ! delta k on discrete wave number
      integer,save      :: NMAX
      real   ,save      :: LFREC,Qq
      complex*16, dimension(:), allocatable :: expK
      complex*16, dimension(:,:), allocatable :: Vout
      end module waveNumVars
      
      module refSolMatrixVars
      !reference solution variables:
      complex*16 :: cOME  
      complex*16, save, allocatable :: B(:)
      complex*16, dimension(:,:,:), allocatable :: Ak
      complex*16, dimension(:), allocatable :: this_B                             
      integer, dimension(:), allocatable :: IPIV
      integer :: info    
      end module refSolMatrixVars
      
      module waveVars
!     save
      real, save :: Escala
      real, save :: Theta !grados
      real, save :: Dt  !segundos
      real, save :: t0
      complex*16, allocatable :: auxUo(:)
      complex*16, dimension(:), allocatable :: Uo
      real, save :: Ts 
      real, save :: Tp 
      real, save :: Zstart
      end module waveVars
            
      module GeometryVars
      ! polynomial fit parameters:
!     use gloVars, only: dp
      integer,  parameter :: degree = 4
      real,save,  dimension(degree+1) :: surf_poly 
      
      !integer, save :: numberOffixedBoundaryPoints
      integer, save :: nXI !number of subdivided segment nodes
      real,save, dimension(:), allocatable :: x,y !subdivided nodes
      
      !                                             ,--Original nodes
      real,save, dimension(:,:,:), allocatable :: Nodes_xz_n
      real,save, dimension(:,:), allocatable :: Xcoord 
            
      ! length and normals at midpoint of BigSegments:
      real, allocatable :: midPoint(:,:)
      real,save, allocatable :: normXI(:,:),lengthXI(:),layerXI(:)
      end module GeometryVars
      
      module resultVars
!     use gloVars, only: dp
       type Pointe
        real :: X,Z
       end type Pointe
       
       type MecaElem
        type (Pointe) :: center
        complex*16, dimension(5) :: Rw !u1,u3,s33,s31,s11
       end type MecaElem
       
       type tipo_qlmk
      !            (imec) l ---, ,--- m (1 = horz, 2= vert)
      !           1 / iGQ ---, | | ,--- k
      !                      | | | | 
       complex*16, dimension(:,:,:,:), allocatable :: qlmk
       logical :: shouldUseQuadrature
       end type tipo_qlmk
       
       type Punto
        !                 ,--- 1 for x, 2 for z
        real    :: center(2)
        real    :: bord_A(2),bord_B(2)
        real    :: normal(2)
        real    :: length
        integer :: layer
        logical :: isBoundary 
        logical :: isOnInterface
        logical :: guardarFK
        logical :: guardarMovieSiblings
        
      !                       ,--- f: 1...nfrec+1
      !                       | ,--- k: 1...NMAX+1 / 2*NMAX
      !                       | | ,--- iMec: 1:5
      !                       | | | 
        complex*16, dimension(:,:,:), allocatable :: FK
        complex*16, dimension  (:,:), allocatable :: FKh,FKv
        
      ! para las funciones de Green:
      !                            ,--- xi
      !                            |  
        type(tipo_qlmk), dimension(:), allocatable     :: GT_k
        
      !               (imec) l ---, ,--- m (1 = horz, 2= vert)
      !                  iGQ ---, | |   ,--- k
      !                 xi ---, | | |   | 
      !                       | | | |   | 
!       complex*16, dimension(:,:,:,:,  :), allocatable :: GT_gq_k
        complex*16, dimension(:,:,:,:)    , allocatable :: GT_gq
!       complex*16, dimension(:,:,:,:,:,:,:), allocatable :: GT_gq_mov_k
        complex*16, dimension(:,:,:,:,:)  , allocatable :: GT_gq_mov
      !                               |
      !                               '--- iMov
        
        
        ! espectro campo total inquirePoints : 
      !                         ,--- f: 1...nfrec+1
      !                         | ,--- iMec: 1:5
        complex*16, dimension  (:,:), allocatable :: W 
        ! espectro campo total moviePoints : 
        complex*16, dimension(:,:,:), allocatable :: WmovieSiblings
      !                       |
      !                       `--- sibling index
      
      !               ,--- xXx (indice punto integracion Gaussiana)
      !               | ,--- (1,2) -> (y,z)
      !               | |
      real, dimension(:,:), allocatable :: Gq_xXx_coords
      real, dimension(:), allocatable :: Gq_xXx_C
      
       end type Punto   
      ! bondary elements:     ,--- POINT index / x (receptor)
      type (Punto), dimension(:), allocatable, save :: allpoints
      type (Punto), dimension(:), allocatable, save :: inqPoints
      type (Punto), dimension(:), allocatable, save :: moviePoints
      type (Punto), dimension(:), allocatable, save :: BouPoints !xi
            
      integer, save :: nIpts, nMpts, NxXxMpts, nBpts, nPts,&
                       iPtini,iPtfin,mPtini,mPtfin,bPtini,bPtfin
      
      !                     ,-,--- nBpts x nBpts (considering 1 direction)
      complex*16, dimension(:,:), allocatable :: ibemMat
      complex*16, dimension(:), allocatable :: trac0vec 
!     complex*16, dimension(:,:,:), allocatable :: ibemPHI
      !                       | `--- m
      !                       `--- i
      integer, dimension(:), allocatable :: IPIVbem
      
      !                                          ,--- ix
      !                                          | ,--- iz
      !                                          | | ,--- iMec
      !                                          | | | ,--- foto{2*Nfrec}
      complex*16, allocatable, save :: Sismogram(:,:,:,:)
      complex*16, dimension(:), allocatable :: auxKvect
!     real*8,dimension(NFREC+1) :: Hf
!     real*8,dimension(NMAX*2) :: Hk
      real*8,dimension(:), allocatable :: Hf,Hk !filtros 
      end module resultVars
              
      module meshVars
!     use gloVars, only: dp
      ! puntos en donde obscultar
      integer, save :: npixX,npixZ
      integer, save :: MeshDXmultiplo, DXmultiploMaxX
      real, save :: MeshDZ,MeshDX,DeltaX
      real, save :: MeshMaxX !from leftmost or rightmost point
      real, save :: MeshMaxZ  
      real, save  :: boundingXn,boundingZn
      real, save  :: boundingXp,boundingZp

      integer, save :: nx,nz
      real, dimension(5,2), save :: colorBounds
      end module meshVars
                  
      module wavelets
      contains 
      
      !  The Ricker wavelet on the time domain saved on   Uo
      subroutine ricker(frec2,outpf)
      use gloVars, only : verbose,PI
      use waveVars, only : Uo,Ts,Tp,Escala,Dt
      implicit none
      integer, intent(in) :: outpf
      integer, intent(in) :: frec2
      integer :: i
      real*8 :: A
      
      ! ajustar el número de puntos a una pontencia entera de 2
!     M = int(log10(real(NUM)) / log10(2.)+0.99 )
!     M = 2**M
      allocate(Uo(frec2))
      Uo=cmplx(0,0,8)
      do i = 1,frec2
        A = pi*(Dt*real(i-1,8)-Ts) / Tp
        Uo(i) = cmplx((A*A-0.5)* exp(- A * A)*Escala,0,8) 
      end do
      if (verbose >= 1) then
       write(outpf,'(A,I4,A,F5.3)') ' A ricker wavelet over a ',frec2, & 
       ' points discrete time signal with Dt = ', Dt
      end if
       
      end subroutine ricker
      
      SUBROUTINE FORK(LX,CX,SIGNI,verbose,outpf)
      implicit none
      integer, intent(in) :: outpf
      integer, intent(in) :: LX,SIGNI,verbose
      COMPLEX*16 :: CARG,CW,CTEMP 
      complex*16,intent(inout) :: CX(LX)
      real*8, parameter :: pi = 4.*ATAN(1.)
      real*8 :: SC
      integer :: i,j,m,istep,l
      if (verbose >= 3) then
        write(outpf,'(a,I4,a)')'FFT on ',LX,' length vector'
      end if
      J=1
      SC=DSQRT(real(1.0,8)/real(LX,8))
      DO 30 I=1,LX
      IF(I > J)GO TO 10
      CTEMP=CX(J)*cmplx(SC,0.0,8)
      CX(J)=CX(I)*cmplx(SC,0.0,8)
      CX(I)=CTEMP
   10 M=LX/2
   20 IF(J <= M)GO TO 30
      J=J-M
      M=M/2
      IF(M >= 1)GO TO 20
   30 J=J+M
      L=1
   40 ISTEP=2*L
      DO 50 M=1,L
      CARG=cmplx(0.0,(pi*real(SIGNI*(M-1)))/real(L),8)  
      CW=EXP(CARG)
      DO 50 I=M,LX,ISTEP
      CTEMP=CW*CX(I+L)
      CX(I+L)=CX(I)-CTEMP
   50 CX(I)=CX(I)+CTEMP
      L=ISTEP
      IF(L < LX)GO TO 40
      
       
      
      RETURN
      END subroutine fork

      end module
      
      module debugStuff
      contains
      
      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(A,E15.5,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E15.5,A)',advance='no') AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      
      subroutine showMNmatrixZabs(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(E15.5,2x)',advance='no') ABS(MAT(i,j))
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      
      !
      subroutine showMNmatrixR(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n ,outpf
      real, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(F10.5,2x)',advance='no') MAT(i,j)
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      !
      subroutine showMNmatrixI(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n, outpf
      integer, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(I10,2x)',advance='no') MAT(i,j)
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      
      end module
      
      !another set of functions in case needed
      module fitting
      contains
      
      function splitatY(surf_poly,degree,Y,aIN,bIN)
      ! there is only one intersection.
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      real, intent(in) :: Y,aIN,bIN
      real, allocatable, dimension(:) :: surf_poly0
      integer, intent(in) :: degree
      integer :: i
      real :: a,b,splitatY,af,bf,temp,tmp2,tmp, & 
                 cx,d,cf,errorTol,s,sf
      logical :: mflag = .true.
      errorTol = real(0.000001,8)
      
      allocate(surf_poly0(degree+1))
      a = aIN
      b = bIN
      
!     write(outpf,*)"orig = ",surf_poly
      surf_poly0 = surf_poly
!     write(outpf,*)"copiado = ",surf_poly0
      surf_poly0(1) = surf_poly0(1) - Y
!     write(outpf,*)"movido = ",surf_poly0
      
      !encontramos el cero de surf_poly0 entre a y b
      af = polyVal(surf_poly0,degree,a)
      bf = polyVal(surf_poly0,degree,b)
!     write(outpf,*)"f(",a,")=",af
!     write(outpf,*)"f(",b,")=",bf
      
            ! usamos el método de Brent.
      if (af*bf >= real(0,8)) then
        if (af < bf ) then
          splitatY = af
        else
          splitatY = bf
        end if
        return
      else
        
        if (abs(af) < abs(bf)) then
          temp = b
          b = a
          a = temp
          temp = bf
          bf = af
          af = temp
        end if
        cx = a
        cf = af
        mflag = .true.
        i = 0
!       write(outpf,*)"a:f(",a,")=",af
!       write(outpf,*)"b:f(",b,")=",bf
        do while((.not.(bf == real(0,8) )) .and. (abs(a-b) > errorTol ))
!          print*,"go_",i
          if ((af /= cf) .and. (bf /= cf)) then
          ! interpolación cuadrática inversa
            s = a * bf * cf / (af-bf) / (af-cf) + & 
            b*af*cf/(bf-af)/(bf-cf)+cx*af*bf/(cf-af)/(cf-bf)
          else
          ! regla de la secante
            s = b - bf * (b-a)/(bf-af)
          end if
          tmp2 = (real(3,8)*a + b)/real(4,8)
          if ( (.not.(((s > tmp2) .and. (s < b)) .or. & 
          ((s < tmp2) .and. (s > b))) ) .or. &
          (mflag .and. ((abs(s-b)) .ge. (abs(b-cx)/real(2,8) ))) .or. &
          ((.not. (mflag)) .and. ((abs(s-b)) .ge. (abs(cx-d)/real(2,8) ))) ) then
            s = (a+b) /real(2,8)
            mflag = .true.
          else
            if ((mflag .and. (abs(b-cx)< errorTol)) .or. &
            ((.not. (mflag)) .and. (abs(cx-d) < errorTol))) then
              s = (a+b) / real(2,8)
              mflag = .true.
            else
              mflag = .false.
            end if
          end if
           sf = polyVal(surf_poly0,degree,s)
           d = cx
           cx = b
           cf = bf
           if (af * sf < real(0,8)) then
             b = s
             bf = sf
           else
             a = s
             af = sf
           end if
           if (abs(af) < abs(bf)) then
             tmp = a
             a = b
             b = tmp
             tmp = af
             af = bf
             bf = tmp
           end if
           i = i + 1
           if (i> 1000) then
             !error
             b = real(0.123456789,8)
             exit
           end if
        end do
      splitatY = b
      end if
      
      
!     splitatY = polyVal(surf_poly,degree,0.)
! print*,"val at ",splitatY," is =",polyVal(surf_poly,degree,splitatY)
! print*,"where,",splitatY,"is =",polyVal(surf_poly0,degree,splitatY),"= 0"
      
      end function 
      
      
      ! function evaluation
      function polyVal(surf_poly,degree,X)
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      integer, intent(in) :: degree 
      real, intent(in) :: X
      real :: polyVal
      integer :: i
      
      !surfo_poly are the polynomial coefficients: A0 A1 A2...An
      polyVal = surf_poly(1) !A0
      DO i = 1,degree
      polyVal = polyVal + X**i * surf_poly(i+1)
      end do
      
      end function
      
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
      ! http://rosettacode.org/wiki/Polynomial_regression#Fortran
      function polyfit(vx,vy,Ln,d,verbose,outpf)
      implicit none
      integer, intent(in)             :: verbose,outpf
      integer, intent(in)                   :: Ln, d
      integer, parameter             :: dp = selected_real_kind(15, 307)
      real, dimension(d+1)              :: polyfit
      real,     dimension(:), intent(in)    :: vx, vy
!     real(dp), dimension(Ln)               :: vx, vy
      real, dimension(:,:), allocatable :: X
      real, dimension(:,:), allocatable :: XT
      real, dimension(:,:), allocatable :: XTX
      integer :: i, j
      integer     :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      real, dimension(:), allocatable :: work
      
      n = d+1
      lda = n
      lwork = n
 
      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, Ln))
      allocate(X(Ln, n))
      allocate(XTX(n, n))
      
      if (verbose >= 2) then
        write(outpf,*)"fitting curve.."
      end if !
      if (verbose >= 3) then
       write(outpf,'(a)')"begin curve fit with:"
       write(outpf,'(a,I4)')"vx: of size:",size(vx)
       write (outpf, '(F9.4)') vx
       write(outpf,'(a,I4)')"vy: of size:",size(vy)
       write (outpf, '(F9.4)') vy
      end if
      
      ! prepare the matrix
      do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
      end do
      if (verbose >= 3) then
       write(outpf,'(a,I4,a,I4,a)') & 
                              "X: size, (",size(X,1),',',size(X,2),')'
       write(outpf,'(10F9.4)') ((X(i,j),j=1,size(X,1)),i=1,size(X,2)) 
      end if
      
      XT  = transpose(X)
      XTX = matmul(XT, X)
      if (verbose >= 3) then
       write(outpf,'(a,I4,a,I4,a)') & 
                         "XTX: size, (",size(XTX,1),',',size(XTX,2),')'
       write(outpf,*)XTX
      end if
      
      ! calls to LAPACK subs DGETRF and DGETRI
      ! factorizacion LU
      call SGETRF(n, n, XTX, & 
                  lda, ipiv, info)
      if (verbose >= 3) then
       write(outpf,'(a)')"DGETRF (XTX):"
       write (outpf, '(E12.3)') XTX
      end if
      !
      if ( info /= 0 ) then
       write(outpf,*) "problem DGETRF =",info
       stop 1
      end if
      ! inversa de la matriz 
      call SGETRI(n, XTX, lda, ipiv, work, lwork, info)
      if (verbose >= 3) then
       write(outpf,'(a)')"DGETRI (XTX):"
       write (outpf, '(F9.4)') XTX
      end if
      
      if ( info /= 0 ) then
       write(outpf,'(a)') "problem DGETRI =",info
       stop 1
      end if
 
      polyfit = matmul( matmul(XTX, XT), vy)
      
      
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      if (verbose >= 2) then
       write(outpf,'(a)')'polyfit='
       write (outpf, '(E12.3)') polyfit
      end if
      end function
      
      !normal vectors to surf_poly at the nXI points (x,y) 
      function normalto(x,y,nXI,surf_poly,degree,verbose,outpf)
      implicit none
      integer,intent(in)      :: verbose,outpf
      integer, intent(in)     :: degree 
!     integer, parameter      :: dp = selected_real_kind(15, 307)
      real, intent(in), dimension(degree+1) :: surf_poly !surface 
      integer :: nXI !number of points
      real, intent(in), dimension(nXI) :: x,y !points coordinates
      
      integer :: i,j
      real, allocatable :: normalto(:,:) !function result value
      real, parameter      :: a = real(1,8) ! normal poing up parameter
      real, parameter      :: tolerance = real(0.001,8) !to tell if its vertical
      real*8, dimension(degree) :: fprime !surf_poly derivative coeficients
      real :: fprimeX, x0, mag 
      
      ALLOCATE (normalto(nXI,2))
      if(verbose >= 2)then
       write(outpf,'(a)',advance='no') "  ...getting the normal vectors"
      end if
      ! the normal line to curve f(x) = An x^n + ... + A1 x + A0 
      ! is y = y1 - 1 / (f'(x1)) (x - x1) 
      ! if we want the normal vectors pointing up (bigger 'y')
      ! we force   y = y1 + a   and we find x in terms of x1
      
      ! f'(x) coefficients
      do i = degree,1,-1
!      write(outpf,*)'i=',i
       fprime(i)= real(i,8) * surf_poly(i+1) 
      end do
      if (verbose >= 2 ) then
       write(outpf,'(a)')'surf_poly=A0,A1,...,AN '
       write(outpf,'(F9.4)') surf_poly
       write(outpf,'(a)')'f_prime_(x) A0 A1 ... An '
       write(outpf,'(F9.4/)') fprime
      end if
      
      do i = 1,nXI !for every point
      !the derivative at (xi)
       fprimeX = real(0,8)
       do j = size(fprime),2,-1
         fprimeX = fprimeX + fprime(j) * x(i)**(j-1)
       end do
       fprimeX = fprimeX + fprime(1)
       x0 = x(i) - a * fprimeX
       
       !normalizar y poner como vector
       mag = sqrt((x(i)-x0)**2 + a**2)
       
       normalto(i,1)= (x0-x(i))/mag
       normalto(i,2)= a/mag
       
       if (verbose >= 3) then
         write(outpf,'(A,f7.2,A,f7.2,A,/A,f6.1,A,f6.1,/A,f6.2,A,f6.2)') &
         "(",x(i),",",y(i),")","x0= ",x0,"  y0= ",y(i)+a, &
         "nx=",normalto(i,1),"  ny=",normalto(i,2)
       end if
       
       
       
      end do
      
      if(verbose >= 2)then
       write(outpf,'(a)',advance='yes') " ... done"
      end if
      
      ! if by incrementing  y = y1 + a the value becomes infinite, then
      ! it is a vertical surface.                                     TODO
      
      
      end function
      
      end module
      module ploteo10pesos
      contains
      
      subroutine plotSpectrum(y_in,Df,full_n,n,titleN,xAx,yAx,logflag,W,H)
      ! (Uo,DFREC,size(Uo),size(Uo)/2.0,titleN,xAx,yAx,logflag,1200,800)
      use DISLIN
      implicit none
      real, intent(in)                              :: Df
      integer, intent(in)                           :: full_n,n,H,W
      character(LEN=9), intent(in)                 :: xAx,yAx
      character(LEN=9)                             :: logflag
      character(LEN=100)                            :: titleN
      COMPLEX*16, DIMENSION(full_n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      
      real, dimension(n) :: x
      real maxY,minY,maxYc,minYc,xstep,ystep
      integer :: i
      character(LEN=100) :: dumb
      CHARACTER(LEN=6)  :: CBUF
      character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      DO i = 1,n
        x(i) = Df*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i))) 
!       !write(6,*) x(i), y(i)
      END DO
      
!     
      minY=MINVAL(real(y(:)),1)
!     !write(6,*)"MinReal Val= ",minY
      maxY=MAXVAL(real(y(:)),1)
!     !write(6,*)"MaxReal Val= ",maxY
      minYc=MINVAL(aimag(y(:)),1)
!     !write(6,*)"MinComplex Val= ",minYc
      maxYc=MAXVAL(aimag(y(:)),1)
!     !write(6,*)"MaxComplex Val= ",maxYc
      minY =MIN(minY,minYc)
      maxY =MAX(maxY,maxYc)
!     !write(6,*)"MinY Val= ",minY
!     !write(6,*)"MaxY Val= ",maxY
      
      logflag = trim(adjustl(logflag))
      if (trim(logflag) == 'logx') then
       logflag = 'X'
       ! los ceros no son muy populares:
       x(1)=x(2)/2.
      elseif (trim(logflag) == 'logy') then
       logflag = 'Y'
!     ! minY =MIN(minY,minYc,-0.1)
!     ! maxY =MAX(maxY,maxYc, 0.1)
      elseif (trim(logflag) == 'loglog') then
       logflag = 'XY'
       ! los ceros no son muy populares:
       x(1)=x(2)/2.
!     ! minY =MIN(minY,minYc,-0.1)
!     ! maxY =MAX(maxY,maxYc, 0.1)
      else
       logflag = 'none'
      end if
      
      
      
      ! Dislin plotting routines 
      CALL METAFL('PDF') !define intended display  XWIN,PS,EPS,PDF
!     !write(titleN,'(a,a)') trim(titleN),'.eps' 
!     !titleN = trim(titleN)
!     !titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
!     !CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL HWFONT()
      CALL axspos (int(370,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W,4) ,int(H,4)) !size of the axis system.
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(5,4) ,'XY') 
!     !call titlin ( titleN , 1 )
      xstep = x(n) /6. ! increment for labels 
      ystep = (maxY-minY)/6.0
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      
      !
      if (logflag == 'Y') then
        CALL AXSSCL('LOG',trim(logflag))
!       !print*,x(n)
!       !print*,log10(x(n))
        call graf(real(x(1),4),real(x(n),4) &
             ,real(x(1),4) ,real(xstep,4) &
             ,real(-7.0,4),real(log10(maxY),4) &
             ,real(-7.0,4),real(1.0,4))  
             
      elseif (logflag == 'X' .or. logflag == 'XY') then
        CALL AXSSCL('LOG',trim(logflag))
        call graf(real(-2.0,4),real(log10(x(n)),4) &
             ,real(-2.0,4) ,real(1.0,4) &
             ,real(minY,4),real(maxY,4) &
             ,real(minY,4),real(ystep,4))
      else
        call graf(real(0.0,4),real(x(n),4), & 
                  real(0.0,4),real(max(1.0,xstep),4), &
                  real(minY,4),real(maxY,4), & 
                  real(minY,4),real(ystep,4))
      end if
      
      
      call color ('RED')
      call curve(real(abs(x),4) ,real(y,4) ,int(n,4))
      call color('BLUE')
      call curve(real(abs(x),4), real(aimag(y),4), int(n,4))
      call color ('FORE') 
      call curve(real(abs(x),4), & 
                    real(sqrt(real(y)**2.0 +aimag(y)**2.0),4), int(n,4))
!     
!     
!     
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
!     
!     
      call legini(CBUF,int(3,4),int(20,4))
!     !nx = nxposn(x(n) + x(n) / 20.)
!     !ny = nyposn(minY + (maxY-minY)*0.7)
!     !print*,nx
!     !print*,ny
      call legpos(int(1600,4),int(320,4))
      write(dumb,'(a)') 'Re(z) '
!     !print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      write(dumb,'(a)') 'Im(z) '
!     !print*,dumb
      call leglin(CBUF,dumb,int(2,4))
      write(dumb,'(a)') 'Abs(z)'
!     !print*,dumb
      call leglin(CBUF,dumb,int(3,4))
      
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(3,4))
!     
      call errmod ("all", "off") !suppress dislin info
      call disfin()
      
!     !print*,'plotted ',trim(titleN)
!     !print*,''
      
      end subroutine plotSpectrum
      
      
      subroutine plotXYcomp(y_in,Dt,n,titleN,xAx,yAx,W,H)
      ! (Uo,Dt,size(Uo),'FIGURE_NAME.pdf','time[sec]','amplitude',1200,800) 
      USE DISLIN
      implicit none
      real, intent(in)                              :: Dt
      integer, intent(in)                           :: n,H,W
      character(LEN=9), intent(in)                 :: xAx,yAx
      character(LEN=100)                            :: titleN
      COMPLEX*16, DIMENSION(n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      
      real, dimension(n) :: x
      real maxY,minY,maxYc,minYc,xstep,ystep
      integer :: i
!     integer :: Modo,nx,ny
      character(LEN=100) :: dumb
      CHARACTER(LEN=30) :: CBUF
      character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      
!     allocate(x(n))
!     print*,size(y)
      DO i = 1,n
        x(i) = Dt*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i))) 
!       write(6,*) x(i), y(i)
      END DO
      
      minY=MINVAL(real(y(:)),1)
!     write(6,*)"MinReal Val= ",minY
      maxY=MAXVAL(real(y(:)),1)
!     write(6,*)"MaxReal Val= ",maxY
      minYc=MINVAL(aimag(y(:)),1)
!     write(6,*)"MinComplex Val= ",minYc
      maxYc=MAXVAL(aimag(y(:)),1)
!     write(6,*)"MaxComplex Val= ",maxYc
      minY =MIN(minY,minYc)
      maxY =MAX(maxY,maxYc)
      
      
!     print*,"plotting"
! Dislin plotting routines ...
      CALL METAFL('PDF') !define intended display XWIN,PS,EPS,PDF
!     write(titleN,'(a,a)') trim(titleN),'.eps' 
!     titleN = trim(titleN)
!     titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
!     print*,"file: ",trim(adjustl(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
!     print*,"disini"
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL HWFONT()
      CALL axspos (int(370,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W,4) ,int(H,4)) !size of the axis system.
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(5,4) ,'XY') 
!     call titlin ( titleN , 1 )
      xstep = x(n)/6.0 ! incremen for labels
      ystep = (maxY-minY)/6.0
      
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0 !solo 3 etiquetas
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      call graf(real(x(1),4), & 
                real(x(n),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
      
!     call title 
      call color ('RED') 
      call curve(real(x,4) ,real(y,4) ,int(n,4))
      call color('BLUE')
      call curve(real(x,4), real(aimag(y),4), int(n,4))
      call color ('FORE') 
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
      call legini(CBUF,int(2,4),int(20,4))
 !     nx = nxposn(x(n)*n + x(n)*n / 20.)
 !     ny = nyposn(minY + (maxY-minY)*0.7)
 !     print*,nx
 !     print*,ny
      call legpos(int(1600,4),int(320,4))
      write(dumb,'(a)') 'Re(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      write(dumb,'(a)') 'Im(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(2,4))
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(2,4))

      call errmod ("protocol", "off") !suppress dislin info
      call disfin()      
      
!     print*,'plotted ',trim(titleN)
      end subroutine plotXYcomp
      
      subroutine plotXYabs(y_in,Dt,n,titleN,xAx,yAx,W,H)
      ! (Uo,Dt,size(Uo),'FIGURE_NAME.pdf','time[sec]','amplitude',1200,800) 
      USE DISLIN
      implicit none
      real, intent(in)                              :: Dt
      integer, intent(in)                           :: n,H,W
      character(LEN=9), intent(in)                 :: xAx,yAx
      character(LEN=100)                            :: titleN
      COMPLEX*16, DIMENSION(n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      
      real, dimension(n) :: x
      real maxY,minY,xstep,ystep
      integer :: i
!     integer :: Modo,nx,ny
      character(LEN=100) :: dumb
      CHARACTER(LEN=30) :: CBUF
      character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      
!     allocate(x(n))
!     print*,size(y)
      DO i = 1,n
        x(i) = Dt*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i))) 
!       write(6,*) x(i), y(i)
      END DO
      
!     minY=MINVAL(real(y(:)),1)
!     write(6,*)"MinReal Val= ",minY
!     maxY=MAXVAL(real(y(:)),1)
!     write(6,*)"MaxReal Val= ",maxY
!     minYc=MINVAL(aimag(y(:)),1)
!     write(6,*)"MinComplex Val= ",minYc
!     maxYc=MAXVAL(aimag(y(:)),1)
!     write(6,*)"MaxComplex Val= ",maxYc
      minY =MINval(abs(y(:)))
      maxY =MAXval(abs(y(:)))
      
      
!     print*,"plotting"
! Dislin plotting routines ...
      CALL METAFL('PDF') !define intended display XWIN,PS,EPS,PDF
!     write(titleN,'(a,a)') trim(titleN),'.eps' 
!     titleN = trim(titleN)
!     titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
!     print*,"file: ",trim(adjustl(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
!     print*,"disini"
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL HWFONT()
      CALL axspos (int(370,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W,4) ,int(H,4)) !size of the axis system.
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(5,4) ,'XY') 
!     call titlin ( titleN , 1 )
      xstep = x(n)/6.0 ! incremen for labels
      ystep = (maxY-minY)/6.0
      
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0 !solo 3 etiquetas
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      call graf(real(x(1),4), & 
                real(x(n),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
      
!     call title 
      call color ('RED') 
      call curve(real(x,4) ,real(abs(y),4) ,int(n,4))
      call color ('FORE') 
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
      call legini(CBUF,int(1,4),int(20,4))
 !     nx = nxposn(x(n)*n + x(n)*n / 20.)
 !     ny = nyposn(minY + (maxY-minY)*0.7)
 !     print*,nx
 !     print*,ny
      call legpos(int(1600,4),int(320,4))
      write(dumb,'(a)') 'abs(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(1,4))

      call errmod ("protocol", "off") !suppress dislin info
      call disfin()      
      
!     print*,'plotted ',trim(titleN)
      end subroutine plotXYabs
      
      subroutine plotFK(thisFK,form,coords,tt,xAx,yAx,outpf)
      !                          |
      !                          `--- 1 : real
      !                          `--- 2 : imag
      !                          `--- 3 : abs with enhanced sharpness 
      use DISLIN
      use waveNumVars, only : NFREC,NMAX,DFREC,DK
      use glovars
      implicit none
      integer ,intent(in) :: form,outpf
      
      !                     ,--- 0...frec
      !                     | ,--- -k...+k
      !                     | | ,--- iMec: 1:5
      !                     | | | 
      complex*16, dimension(NFREC+1,2*NMAX,imecMax),intent(in) :: thisFK
      real, dimension(2), intent(in) :: coords
      character(LEN=10),intent(in) :: tt
      character(LEN=9), intent(in) :: xAx,yAx
      
      real, dimension(NMAX)             :: vVert
      real, dimension(NFREC)            :: vHorz
      real, dimension(NFREC,NMAX)       :: M
      character(LEN=10), dimension(5)   :: nombre
      character(LEN=100)                :: titulo
      integer :: i,ik,iMec,Sf
!     real :: k
      real :: minX,minY,maxX,maxY,xstep,ystep,miV,maV,Vstep!,x_i,z_i
      real, dimension(41)   :: ZLVRAY
      real, parameter :: p = 12. ! sharpness parameter
      
      
      if(verbose>=1)write(outpf,'(a,a,a)') "will plot ",trim(tt),"..."
            
      M = 0
      nombre(1)= '_w.png'
      nombre(2)= '_u.png'
      nombre(3)= '_s33.png'
      nombre(4)= '_s31.png'
      nombre(5)= '_s11.png'
      
      do i=1,NFREC
         vHorz(i) = 0 + (i-1) * DFREC
      end do
      !
      do ik=1,nmax
         vVert(ik) = 0 + (ik-1) * DK
      end do
      
      minY = 0.
      maxY = vVert(size(vVert))
      ystep = maxY / 10.
      
      minX = vHorz(1)
      maxX = vHorz(size(vHorz))
      xstep = maxX / 10.
      
!     nIP = size(thisFK,4)
!     if(verbose>=1)write(outpf,'(a,I0)')"number of points: ",nIP
      
!     do iP = 1,nIP 
!        x_i=coords(iP,1)
!        z_i=coords(iP,2)
      do iMec = 1,imecMax!min(2,imecMax) !podrian ser más pero meh      
      write(titulo,'(a,a,I0,a,I0,a,I0,a,I0,a,a)') trim(tt),'[', &
      int(coords(1)),'.',abs(int((coords(1)-int(coords(1)))*10)),';', & 
      int(coords(2)),'.',abs(int((coords(2)-int(coords(2)))*10)),']', & 
      trim(nombre(iMec))
      
      if (form .eq. 1) then
       M = real(thisFK(1:NFREC,1:NMAX,iMec))
      elseif (form .eq. 2) then
       M = aimag(thisFK(1:NFREC,1:NMAX,iMec))
      elseif (form .eq. 3) then
       M = log(1. + exp(p)*abs(thisFK(1:NFREC,1:NMAX,iMec))) / & 
           log(exp(p)+1.)
       M = M / maxval(M)
      end if
      
      miV = minval(M)
      maV = maxval(M)
      Vstep = (maV-miV)/40.0
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      
      ! shaded contour plot
      CALL METAFL('PNG') !'PDF'
      CALL SETFIL(trim(titulo))
      if(verbose>=2) write(outpf,'(a)')trim(titulo)
      call filmod('DELETE') ! para sobreescribir el archivo
!     CALL PAGE (2000, 1100)
      Sf = 2
      CALL PAGE (int(1800* SF,4) , int(1200* SF,4))
      
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(800,4))
!     call sclmod('FULL')
!     CALL PAGMOD('NONE')
      CALL SCRMOD('REVERS') !fondo blanco
!     CALL SETPAG('DA4P')
      CALL DISINI()
!     CALL COMPLX() !texto complejo
!     call incmrk ( 1 )
!     CALL HWFONT()
      CALL axspos (int(180* SF,4) ,int(1100* SF,4)) !the position of an axis system. Lower left corner
      call axslen (int(1400* SF,4),int(1000* SF,4)) !size of the axis system.
      call labdig(int(2,4),'XY') !number of decimal places for labels
      call labdig(int(1,4),'Z') !number of decimal places for labels
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y')
      
!     CALL LABELS ('EXP ','Z')
      CALL SETVLT ('SPEC')
      
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(minY,4),real(maxY,4),real(minY,4),real(ystep,4), & 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 

      
      CALL CONSHD(real(vHorz,4), int(size(vHorz),4), & 
                  real(vVert,4), int(size(vVert),4), & 
                  real(M,4), real(ZLVRAY,4),int(41,4))
      
      CALL HEIGHT(int(50,4))
      call errmod ("all", "off")
      CALL DISFIN()
      
      end do !iMec
!     end do !iP
      end subroutine plotFK
      
      
      
      subroutine drawGEOM(titleN,whole,outpf)
      use DISLIN
      use soilVars, only : Z,N
      use GeometryVars, only: nXI,Xcoord
      use glovars
      implicit none
      logical :: whole
      character(LEN=100) :: titleN!,dumbTxt
      real :: maxY,minY,maxX,minX,xstep,zstep
!     real :: xi(nXI),yi(nXI)
      integer :: J,outpf
      if (verbose >= 2) Write(outpf,'(a)') "Will plot geometry..."
      
      ! plotting boundaries
      minX = MINVAL(Xcoord(:,1))
      minX = MIN(minX,-1.)
!     print*,"minX ",minX
      maxX = MAXVAL(Xcoord(:,1))
      maxX = MAX(maxX, 1.)
!     print*,"maxX ",maxX
!     minY = MINVAL(y)
      minY = MIN(MINVAL(Xcoord(:,2)), 1.)
!     print*,"minY ",minY
      
      if (whole .eqv. .true.)then
        maxY = max(maxval(Xcoord(:,2))+maxval(Xcoord(:,2))/10,Z(N+1))
      else
        maxY = max(maxval(Xcoord(:,2))+maxval(Xcoord(:,2))/10,1.)
      end if
!     maxY = MAXVAL(y)
!     print*,"maxY ",maxY
!     stop 0
!     maxY = MAX(maxY, 1.)
      
      minX = minX*1.1
      maxX = maxX*1.1
      minY = minY-(maxY/10)
!     minY = MIN(minY,Z(N+1)+10.)
      
      ! Dislin plotting routines 
      CALL METAFL('PDF') ! define intended display  XWIN,PS,EPS,PDF
      CALL SETFIL(trim(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(1200+300,4),int(800+300,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      call incmrk (int(1,4))
      CALL HWFONT()
      CALL axspos (int(250,4) ,int(800+150,4)) !the position of an axis system. Lower left corner
      call axslen (int(1200,4) ,int(800,4)) !size of the axis system.
!     call name(trim(xAx),'X') 
!     call name(trim(yAx),'Y') 
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(10,4) ,'XY') 
      !            low X   left Y   upp X   right Y
      call setgrf("NAME", "NAME", "NONE", "LINE")
      CALL SETVLT ('SPEC')
!     call titlin ( titleN , 1 )	
      ! increment for labels 
      xstep = real(   ( (maxX-minX) / 5.0 ))
      zstep = real(   ( (maxY-minY) / 10. ))
!     print*,minX
!     print*,maxX
!     print*,xstep
!     print*,minY
!     print*,maxY
!     print*,(maxY-minY) / 10.
!     print*,zstep
!     stop 0
      ! xini xfin yfirstvalue xstep ymin ymax yminvalueTicks yTicks
      call graf(real(minX,4),real(maxX,4),real(minX,4),& 
                real(xstep,4),real(maxY,4),real(minY,4),& 
                real(maxY,4),real(-zstep,4)) 
      
      do J=1,N
         call color ('FORE')
         call shdpat(int(J,4))
         call rlrec(real(minX,4),real(Z(J),4),& 
                    real(maxX-minX,4),real(Z(J+1)-Z(J),4))
      end do
      
      call color ('BACK')
      call shdpat(int(16,4))
      call rlarea(real(Xcoord(:,1),4),real(Xcoord(:,2),4),int(nXI,4))
      call color ('FORE')
      
      CALL HSYMBL(int(7,4))
      CALL MRKCLR(int(-1,4))
      call marker(int(15,4))
      
      call curve(real(Xcoord(:,1),4),real(Xcoord(:,2),4),int(nXI,4))
      call color ('BLACK')
      call rline(real(Xcoord(1,1),4),real(Xcoord(1,2),4),real(Xcoord(nXI,1),4),real(Xcoord(nXI,2),4))
      
!     call xaxgit() 
      call errmod ("all", "off")
      call disfin()
      
      end subroutine drawGEOM
      
      
      
      end module ploteo10pesos
      PROGRAM reference_solution
      use refSolMatrixVars ! cOME,A,B,XX,IPIV,info
      use gloVars !UI,UR,PI,verbose
      use waveNumVars !NFREC,FREC,DFREC,OME,OMEI,DK,NMAX,LFREC,DK,NMAX,expk
      use soilVars, only : N
      use debugStuff
      use resultVars
      use GeometryVars!, only : numberOffixedBoundaryPoints
      use ploteo10pesos
      use sourceVars
      use wavelets
      use MeshVars
      use Gquadrature, only: Gquad_n
      use waveVars, only : Uo,Escala
      
      implicit none
      ! output direction   6 : screen      101: log file
      integer, parameter :: PrintNum = 6
      !Loop Counters
      integer :: J  ! frequency loop
      integer :: l,m,iP,iPxi,iP_x,i,ik !small loop counter
      integer :: status 
      CHARACTER(len=400) :: path
      character(LEN=100) :: titleN
      character(LEN=10) :: tt
      character(LEN=9)  :: xAx,yAx
      real :: factor,k
      complex*16 :: integralEq16,integralEq16mov,freefTraction
      
      character(10) :: time
      real :: startT,finishT
      call cpu_time(startT)
!     complex*16 :: integralEq14
      ! output direction
      if (PrintNum /= 6) open(PrintNum,FILE= "GlayerOUT.txt")
      call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
      CALL getcwd(path)
      write(PrintNum,'(a,/,a)') 'At:',TRIM(path)
      write(path,'(a,a)') trim(path),"/WorkDir"
      CALL chdir(trim(path))
      CALL getcwd(path)
      write(PrintNum,'(a,/,a)') 'Now at:',TRIM(path) 
      
      CALL chdir("outs",status)
      if (status .eq. 0) call system("rm *.*")
      if (status .eq. 0) call chdir("..")
      
      call getMainInput
      
      call getSoilProps (PrintNum)    
      nIpts=0; nMpts=0; nBpts = 0
      iPtfin = 0
      mPtfin = 0
      allocate(vout(1,2))
      if (getInquirePointsSol) call getInquirePoints(PrintNum)
      call getVideoPoints(PrintNum)
      if (workBoundary) call getTopography(PrintNum)
      NPts = nIpts + nMpts
      allocate (allpoints(Npts))
      allocate (Ak(4*N+2,4*N+2,nmax+1))
      ALLOCATE (B (4*N+2)) ! una sola fuente
      ALLOCATE (IPIV(4*N+2)) ! pivote
      allocate (auxKvect(2*Nmax))
      allocate(expK(2*nmax))
      factor = sqrt(2.*NMAX*2)
      
      if (getInquirePointsSol) then 
         allpoints(iPtini:iPtfin)= inqPoints
         deallocate(inqPoints); end if!
      if (makeVideo) then
         allpoints(mPtini:mPtfin) = moviePoints
         deallocate(moviePoints); end if
      
      ! Espectro - numero de onda
      do iP=1,Npts; if (allpoints(iP)%guardarFK)then
          allocate(allpoints(iP)%FK(NFREC+1,NMAX,5))
      end if; end do!
      
      do iP=iPtini,iPtfin !solo inquirepoints:
         allocate(allpoints(iP)%W(NFREC+1,5)); allpoints(iP)%W = 0
      end do
      
      if (makeVideo) then; do iP=mPtini,mPtfin
         allocate(allpoints(iP)%WmovieSiblings(NxXxMpts,NFREC+1,5))
         allpoints(iP)%WmovieSiblings = 0
      end do; end if!
      
      if (workBoundary) then
      ! G para resolver el campo difractado por topografía
       do iP_x = 1,nPts 
         allocate(allpoints(iP_x)%GT_k(nBpts))
       end do!
       do iP_x = iPtini,iPtfin
         allocate(allpoints(iP_x)%GT_gq(nBpts,Gquad_n,5,2))
       end do
       if (makeVideo) then
         do iP_x = mPtini,mPtfin
           allocate(allpoints(iP_x)%GT_gq_mov(nBpts,Gquad_n,5,2,NxXxMpts))
         end do
       end if
      end if
      
      ! source application point:
      xfsource = allpoints(1)%center(1)
      zfsource = allpoints(1)%center(2)
      efsource = allpoints(1)%layer
      intfsource = allpoints(1)%isOnInterface
      nxfsource = allpoints(1)%normal(1)
      nzfsource = allpoints(1)%normal(2)
      
      if (verbose .ge. 1) then
        write(PrintNum,'(a,F8.2,a,F8.2,a)') & 
                            "source is at (",xfsource,",",zfsource,")"
      end if
      
      call makeTaperFuncs(20,0.8,0.8)
      call expPIK(expK)
      call waveAmplitude(PrintNum)
      DO J=1,NFREC+1
        ! complex frecuency for avoiding poles and provide damping
        FREC=DFREC*real(J-1)
        if (J .eq. 1) FREC = DFREC*0.001 
        OME=2.0*PI*FREC

!       OMEI=0.7*PI/TW
!       OMEI=DFREC / 2.0
        ! Bouchon (2003) OMEI entre -pi/T y -2pi/T ; T= 2pi/DFREC
!       COME=CMPLX(OME, -OMEI*1.0) !periodic sources damping
!       COME=COME*(UR - UI/2.0/Qq) !histeretic damping
        cOME = cmplx(OME, -DFREC / periodicdamper,8) * cmplx(1.0, -1.0/2.0/Qq,8)
        
        if(verbose>=1)then
 write(PrintNum,'(A,I0,A,EN13.2,1x,EN13.2,A,EN11.2,a,EN11.1,a)') &
 'w(',J,')= ',REAL(COME),aimag(COME),'i :: ',FREC,' Hz :: ',OME,' rad/s'
        end if 
        
      ! Subsegment the topography if neccesssary
      if (workBoundary) then 
         call subdivideTopo(J,PrintNum)
         call allocintegPoints(J)
      end if
      
      Ak = 0
      Do ik=1,nmax+1 !WAVE NUMBER LOOP-
         k = real(ik-1,8) * dK
         if (ik .eq. 1) k = dk * 0.001
         call matrixA_borderCond(Ak(:,:,ik),k,cOME,PrintNum)
         call inverseA(Ak(:,:,ik),4*N+2)
      end do ! ik
      
        ! for the discrete wave number
!     LFREC=2000 + L*exp(-PI*(FLOAT(J-1)/NFREC)**2)  !EMPIRICAL (improve)
!     DK=2.0*PI/LFREC
      
!     NMAX=(OME/minval(BETA))/DK+1                   !EMPIRICAL (improve)
!     NMAX=2*NMAX+100                                !EMPIRICAL (improve)
      
      call stratumDiffField (0, &
                             allpoints, npts, 0, &
                             J, cOME, PrintNum, Vout)
      
!     do ip_x = iPtini,iPtfin
!        allpoints(ip_x)%W(J,:) = allpoints(ip_x)%W(J,:) * Uo(J) * Escala
!     end do
!     if (makeVideo) then
!      do ip_x = mPtini,mPtfin
!        allpoints(ip_x)%WmovieSiblings(:,J,:) = allpoints(ip_x)%WmovieSiblings(:,J,:) * Uo(J) * Escala
!      end do
!     end if
      !
      if (workboundary) then
      trac0vec = cmplx(0.0,0.0,8) 
      ! diffracted field by stratification at the boundary
      call stratumDiffField (0, &
                             boupoints, nBpts, 0, &
                             1, cOME, PrintNum, Vout)
!     do ip_x = 1,nbpts
!        boupoints(ip_x)%W(1,:) = boupoints(ip_x)%W(1,:) * Uo(J) * Escala
!     end do
      do iP_x = 1,2*nBpts,2
       do l=0,1 !dirección de tracción en el receptor (mini renglón)
         trac0vec(iP_x + l) = -1.0 * freefTraction(ceiling(iP_x/2.),l)
       end do !l
      end do
      if (verbose .ge. 2) then
        call showMNmatrixZ(2*nBpts,1, trac0vec,"trac0",PrintNum)
      end if
      
      ! llenar la matriz de ibem por columnas
      ibemMat = cmplx(1.0,0.0,8) 
      do iPxi = 1,2*nBpts,2
        vout = cmplx(0.0,0.0,8)
        call stratumDiffField (ceiling(iPxi/2.), &
                             boupoints, nBpts, -1, &
                             J, cOME, PrintNum, Vout)
        ibemMat(:,iPxi:iPxi+1) = Vout
        
        call stratumDiffField (ceiling(iPxi/2.), &
                              allpoints, npts, 0, &
                              J, cOME, PrintNum, Vout)
      end do
      
      if (verbose .ge. 2) then
        call showMNmatrixZ(2*nBpts,2*nBpts,ibemMat,"ibMat",PrintNum)
      end if
      
      ! resolver ibem
      iPIVbem = 2*nBpts
      call zgesv(2*nBpts,1,ibemMat,2*nBpts,IPIVbem,trac0vec,2*nBpts,info)
      if(info .ne. 0) stop "problem wiht ibem system"
      
      if (verbose .ge. 1) then
       if(verbose .ge. 2) call showMNmatrixZabs(2*nBpts,1, trac0vec," phi ",PrintNum)
         CALL chdir("outs")
         write(titleN,'(a,I0,a)') 'phi_',J,'.txt'
         OPEN(351,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         do i = 1,2*nBpts,2
          write(351,'(ES14.5E2,2x,ES14.5E2,2x,ES14.5E2,6x,ES14.5E2,2x,ES14.5E2,2x,ES14.5E2)') & 
                     real(trac0vec(i)),aimag(trac0vec(i)),abs(trac0vec(i)), &
                     real(trac0vec(i+1)),aimag(trac0vec(i+1)),abs(trac0vec(i+1))
         end do
         close (351) 
         write(titleN,'(a,I0,a)')'phi_',J,'h.pdf'
         call plotXYabs(trac0vec(1:2*nBpts:2),1.0,nBpts,titleN, &
         '  phi h  ','abs(phi) ',1200,800)
         write(titleN,'(a,I0,a)')'phi_',J,'v.pdf'
         call plotXYabs(trac0vec(2:2*nBpts:2),1.0,nBpts,titleN, &
         '  phi v  ','abs(phi) ',1200,800)
         write(titleN,'(a,I0,a)')'phi_',J,'abs.pdf'
         call plotXYabs(cmplx(abs(trac0vec(1:2*nBpts:2))+abs(trac0vec(2:2*nBpts:2)),0.0,8),1.0,nBpts,titleN, &
         '  |phi|  ','abs(phi) ',1200,800)
         CALL chdir("..")
      end if!
      if(verbose .ge. 2) write(PrintNum,'(a)') "add diffracted field by topography"
      ! usar coeficientes de fuerza por segmento 
      ! para encontrar campo difractado por topografía
      do iP_x = iPtini,iPtfin !cada receptor X
        do iPxi = 1,2*nBpts,2 !cada fuente virtual (dos direcciones)
          do l=0,1 !direc. desp. en receptor X  (x,z)
            do m=0,1 !direc. func. Green
!       allpoints(iP_x)%W(J,l+1) = allpoints(iP_x)%W(J,l+1) + &
        allpoints(iP_x)%W(J,2-l) = allpoints(iP_x)%W(J,2-l) + &
        integralEq16(iP_x,l,ceiling(iPxi/2.),m) * trac0vec(iPxi+m)
            end do !m
          end do ! l
        end do ! iPxi
      end do !iP_X inqPoints
      
      if (makeVideo) then
      do iP_x = mPtini,mPtfin ! cada nivel de receptores X
       do iP = 1,NxXxMpts ! cada punto de película
        do l=0,1 !direc. desp. en receptor X  (x,z)
          do iPxi = 1,2*nBpts,2 !cada fuente virtual (dos direcciones)
            do m=0,1 !direc. func. Green
!       allpoints(iP_x)%WmovieSiblings(iP,J,l+1) = allpoints(iP_x)%WmovieSiblings(iP,J,l+1) + &
        allpoints(iP_x)%WmovieSiblings(iP,J,2-l) = allpoints(iP_x)%WmovieSiblings(iP,J,2-l) + &
        integralEq16mov(iP,iP_x,l,ceiling(iPxi/2.),m) * trac0vec(iPxi+m)
            end do !m
          end do !iPxi
        end do ! l
       end do ! iP
      end do !iP_x
      end if
      end if ! workboundary
      END DO ! J: frequency loop
      
!     deallocate(this_B)
      deallocate(B);deallocate(IPIV)
      
      if(verbose >= 1) write(PrintNum,'(a)')"response found..."
      
      ! showoff 
      if (plotFKS) then 
      if(verbose >= 1) write(PrintNum,'(a)')"show off"
           write(tt,'(a)')"FK_"
           write(xAx,'(a)')"frec [Hz]"
           write(yAx,'(a)')" K [1/m] "
           CALL chdir("outs")
         do iP = iPtini,iPtfin
           if (allpoints(iP)% guardarFK) then
      call plotFK(allpoints(iP)%FK,3,allpoints(iP)%center,tt,xAx,yAx,PrintNum)
           end if
         end do
           CALL chdir("..")
      end if
      
      ! mostrar sismogramas en los puntos de interés
           CALL chdir("outs")
           CALL getcwd(path)
           write(PrintNum,'(a,/,a)') 'Now at:',TRIM(path)
      do iP = iPtini,iPtfin
        call W_to_t(allpoints(iP)%W,iP,allpoints(iP)%center,PrintNum)
      end do
           CALL chdir("..")
      
      call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time  
      call cpu_time(finishT)
      write(PrintNum,'(a,f10.3)') "Elapsed time [sec] before video = ",finishT-startT
      
      if (makeVideo) call Hollywood(PrintNum)
      
      Write(PrintNum,'(a)') ' done '    
      END program
      
      subroutine getMainInput
      use glovars
      use Gquadrature, only : WLmulti
      implicit none
      logical :: lexist
      
      inquire(file="maininput.txt",exist=lexist)
      if (lexist) then
        OPEN(35,FILE="maininput.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "maininput.txt" on Working directory' 
      end if
      
      READ(35,*)
      READ(35,'(I1)') verbose; print*,"verbose =",verbose
      READ(35,*) 
      READ(35,'(L1)') makevideo; print*,"make a video =",makevideo
      READ(35,*)
      READ(35,'(L1)') workBoundary; print*,"boundary? ",workBoundary
      READ(35,*)
      READ(35,'(L1)') plotFKS; print*,"plotFK?",plotFKS
      READ(35,*)
      READ(35,*) multSubdiv; print*,"division multiple = ", multSubdiv
      READ(35,*)
      READ(35,*)
      READ(35,*) WLmulti; print*,"integ WL multiple = ", WLmulti
      read(35,*)
      read(35,'(L1)') symmetricProblem; print*,"The problem is symmetric? ", symmetricProblem
      read(35,*)
      read(35,*) periodicdamper; print*,"Periodic sources damping factor ", periodicdamper
      close(35)
      
      end subroutine getMainInput
      subroutine getSoilProps (outpf)
      use soilVars
      use waveNumVars
      use waveVars, only : dt,t0
      use gloVars, only : verbose,PI
      implicit none
      integer, intent(in) :: outpf
      logical :: lexist
      real*8 :: H, ALF, BET, RO
      integer :: J
      
      ! vemos si existen los archivos
      
      inquire(file="HVDEPTH.txt",exist=lexist)
      if (lexist) then
        OPEN(7,FILE="HVDEPTH.txt",FORM="FORMATTED")
      else
        write(outpf,'(a)') 'There is a missing input file. '
        stop 'Check "HVDEPTH.txt" on Working directory' 
      end if
      
      READ(7,*)
      READ(7,'(I2)')N   !NUMBER OF LAYERS; HALF-SPACE AT N+1
      READ(7,*)
      if (verbose >= 1) then
       write(outpf,'(a)')' '
       write(outpf,'(I2,A)') N,' layers'
       write(outpf,'(a)')' '
      end if
      ALLOCATE (Z(N+1))
      ALLOCATE (AMU(N+1))
      ALLOCATE (BETA(N+1))
      ALLOCATE (ALFA(N+1))
      ALLOCATE (LAMBDA(N+1))
      ALLOCATE (RHO(N+1))
      
      Z(1)=real(0,8)
      if (verbose >= 1) then
       write(outpf,'(a)')& 
              '        depth       mu         beta         alpha'
      end if
      DO J=1,N
         READ(7,*) H, ALF, BET, RO
         Z(J+1)=Z(J)+real(H,8)
         AMU(J)=RO*BET**2
         BETA(J)=BET
         ALFA(J)=ALF
         RHO(J) = RO
         LAMBDA(J)=RHO(J)*ALFA(J)**2 - real(2,8)*AMU(J)
!         BEALF=SQRT((0.5-ANU)/(1.0-ANU)) !IF POISSON RATIO IS GIVEN
!         ALFA(J)=BET/BEALF
       if (verbose >= 1) then
         write(outpf,'(F7.1,A,F7.1,2x,F8.2,3x,F7.1,6x,F7.1)') & 
         Z(J),' - ',Z(J+1),AMU(J),BETA(J),ALFA(J)
       end if
      END DO
      
      READ(7,*) H, ALF, BET, RO
      AMU(N+1)=RO*BET**2
      BETA(N+1)=BET
      ALFA(N+1)=ALF
      RHO(J) = RO
      LAMBDA(N+1)=RHO(J)*ALFA(N+1)**2 - real(2,8)*AMU(N+1)
      if (verbose >= 1) then
       write(outpf,'(F7.1,A,5x,F8.2,3x,F7.1,6x,F7.1)') &
       Z(N+1),' - inf.',AMU(N+1),BETA(N+1),ALFA(N+1)
       write(outpf,'(a)')' '
      end if
!      BEALF=SQRT((0.5-ANU)/(1.0-ANU))   !IF POISSON RATIO IS GIVEN
!      ALFA(N+1)=BET/BEALF
      
      READ(7,*)
      READ(7,*)DFREC,NFREC,DK,NMAX,Qq,t0
      close(7)
      
      !Max frequency should be adjusted to a 2**N value
         NFREC= int(log10(real(NFREC)) / log10(2.)+0.99 )
         NFREC= 2**NFREC
         
         NMAX = int(log10(real(NMAX)) / log10(2.)+0.99 )
         NMAX = 2**NMAX
      
!        LFREC = 2*PI/DK
         
      Dt = 1.0 / (2.0 * real(NFREC) * DFREC)
   
      if (verbose >= 1) then
       write(outpf,'(a,I0,a,F8.4,a,F8.4,a,/,a,F8.1,/)') & 
           'N. frequencies: ',NFREC,'  @',DFREC,'Hertz :. Fmax = ', & 
           NFREC*DFREC,'Hertz','Atenuation Q = ',Qq
       
       write(outpf,'(a,I0,/)') 'N. wavenumbers: ',NMAX
        
!       write(outpf,'(a,I0,a,E11.2,a,EN10.2,a,/,a,E11.2,a,E11.2,a,/)') & 
!       'N. wavenumbers: ',NMAX,'  @',dK,' :: @L = ',LFREC,'m',&
!       'k=[',-NMAX*dK,'... 0 ...',NMAX*dK,']'
       
      end if
      
      end subroutine getSoilProps
      subroutine getInquirePoints(outpf)
      use resultVars, only : inqPoints, nIpts, iPtini,iPtfin
      use soilVars, only : Z,N
!     use gloVars, only : dp
      implicit none
      integer, intent(in) :: outpf
      integer :: i,e
      logical :: lexist
      real ::  errT = 0.001
      integer :: auxGuardarFK
      
      ! read file
      inquire(file="interestingPoints.txt", exist=lexist)
      if(lexist)then
        open(7,FILE="interestingPoints.txt",FORM="FORMATTED")
      else
        write(outpf,'(a,a)') 'There is a missing input file. ',&
        'Check "interestingPoints.txt" on Working directory' 
        stop 1
      end if
      READ(7,*) !Points of interest for an accelerogram
      READ(7,*) nIpts
      iPtini = 1
      iPtfin = nIpts
!     bouPtStartInd=numberOfnonBoundaryPoints+1
!     print*,'numpts=',numberOfnonBoundaryPoints
      allocate(inqPoints(nIpts))
      inqPoints(:)%isOnInterface = .false.
      inqPoints(:)%isBoundary = .false.
      inqPoints(:)%guardarFK = .false.
      inqPoints(:)%guardarMovieSiblings = .false.
!     do j=1,N+1
!     print*,"z(",j,")=",Z(j)
!     end do
      
      READ(7,*) !  X        Z          nx       nz     guardarFK
      do i=1, nIpts 
!        allocate(inqPoints(i)%FK(NMAX,NFREC+1,imecMax))
!        inqPoints(i)%FK = 0
      
         READ(7,*) inqPoints(i)%center(1), inqPoints(i)%center(2), &
                inqPoints(i)%normal(1), inqPoints(i)%normal(2), auxGuardarFK
!     print*,'[',inquirePoints(i)%center%X,inquirePoints(i)%center%Z,']'
      
      !encontrar el alyer en el que estan o 0 si está sobre la interfaz
!          print*,"z=",nonBoundPoints(1)%center(i,2),"....."
           do e=1,N
!            print*,"at ",real(Z(e+1))
             if(real(inqPoints(i)%center(2)) .lt. Z(e+1)) then
!              print*,"ok", e
               exit
             end if
           end do
!          print*,"e=",e," --- z=",nonBoundPoints(1)%center(i,2)
           inqPoints(i)%layer = e
           if (auxGuardarFK .eq. 1 ) then
              inqPoints(i)%guardarFK = .true.
           end if
           
           do e=1,N+1
             if((Z(e)-errT .lt. real(inqPoints(i)%center(2))) & 
             .and. (real(inqPoints(i)%center(2)) .lt. Z(e)+errT)) then
!               print*,real(nonBoundPoints(1)%center(i,2))," is on the interface"
                inqPoints(i)%isOnInterface = .true.
             end if
           end do
      end do
!     stop 0
      end subroutine getInquirePoints   
      
      
      subroutine getVideoPoints(outpf)
      use resultVars, only : moviePoints, nMpts, & 
                             iPtfin,mPtini,mPtfin
      use meshVars
      use gloVars,only:verbose,PI,makevideo
      use soilVars, only : Z,N
      use waveNumVars, only : DK,NMAX
      use resultvars, only: NxXxMpts
      use wavevars, only: Escala,Ts,Tp
      implicit none
      integer, intent(in)          :: outpf!,nJ
      logical                      :: isOnIF
      integer                      :: plusCero
      integer                      :: iz,iP,e,currentLayer
      real                         :: errT = 0.001
      logical :: lexist
      
      inquire(file="íncidenceParameters.txt",exist=lexist)
      if (lexist) then
        OPEN(7,FILE="íncidenceParameters.txt",FORM="FORMATTED")
      else
        write(outpf,'(a)') 'There is a missing input file. '
        stop 'Check "íncidenceParameters.txt" on Working directory' 
      end if
      
      READ(7,*)
      READ(7,*)
      READ(7,*) npixZ !num de pixeles
      READ(7,*)
      READ(7,*) npixX !num de pixeles
      READ(7,*)
      READ(7,*) MeshMaxZ
      READ(7,*)
      READ(7,*) MeshMaxX
      READ(7,*)
      READ(7,*)
      READ(7,*) Escala
      READ(7,*)
      READ(7,*) Ts
      READ(7,*)
      READ(7,*) Tp
      close(7)
      
      if (makeVideo .eqv. .false.) return 
        
        write(outpf,'(a,F12.2,a)') "Lfrec = 2*pi/DK = ",2*PI/DK, "m"
        DeltaX = 1.0 / (nMax * DK)
        write(outpf,'(a,F9.3)')"Delta X = ",DeltaX
        boundingXp =  1.0 / DK
        boundingXn = - boundingXp
        if (boundingXp .lt. MeshMaxX) stop "You need to reduce dK"
        ! encontrmos un multiplo de deltaX que sea mayor o igual a meshmaxX
        DXmultiploMaxX = int(MeshMaxX/DeltaX) + 1 ! 11
        ! número de puntos de la película
        NxXxMpts = 2* DXmultiploMaxX + 1 !23
        
        if (npixx .lt. NxXxMpts) then!             10
          MeshDXmultiplo = max(1, int(NxXxMpts / npixX)) ! 2
          NxXxMpts = 2* (int(MeshMaxX/(DeltaX* MeshDXmultiplo)) + 1) + 1
        else
          MeshDXmultiplo = 1
        end if
        
        boundingZp = MeshMaxZ !Z(N+1)
        boundingZn = 0. 
        plusCero = 0
        if(boundingZn == 0.0) then
          plusCero =1
        end if
!       nz = int((abs(boundingZn)+abs(boundingZp)) / MeshDZ) + plusCero
        MeshDZ = MeshMaxZ / npixZ
        nz = npixZ + plusCero
        
        if (verbose .ge. 1) then
          Write(outpf,'(a)')"Video bounding box: "
!         write(outpf,'(a,I0,a)') "nx=(",nx,")"
          write(outpf,'(a,F9.1,a,F9.1,a,F10.3,a,I0,a,I0,a)') & 
                        "x :[", - MeshDXmultiplo * DeltaX * (NxXxMpts-1)/2," ... ", &
                         MeshDXmultiplo * DeltaX * (NxXxMpts-1)/2,"] :: ", DeltaX,"m @[", &
                         MeshDXmultiplo,"]dk = ", NxXxMpts," puntos"
          write(outpf,'(a,F9.1,a,F9.1,a,F10.3,a,I0,a)') & 
                     "z :[",boundingZn," ... ",boundingZp, & 
                     "] :: ",MeshDZ,"m        = ",nz," puntos"
        end if
        
        nMpts = nz
        mPtini = iPtfin + 1
        mPtfin = mPtini + nMpts - 1
        allocate(moviePoints(nMpts))
        moviePoints(:)%isBoundary = .false.
        moviePoints(:)%isOnInterface = .false.
        moviePoints(:)%guardarFK = .false.
        moviePoints(:)%guardarMovieSiblings = .true.
        
        if(verbose>=2)Write(outpf,'(a,I0)') & 
        "Number of verical movie pixels: ", nMpts
        iP = 1
        do iz = 1,nz
         ! el layer para este z
         isOnIF = .false.
         do e=1,N
            if(real(boundingZn + (iz-1) * MeshDZ) .lt. Z(e+1))then
               exit
            end if
         end do
         currentLayer = e
         do e=1,N+1
            if((Z(e)-errT .lt. real(boundingZn + (iz-1) * MeshDZ)) & 
         .and. (real(boundingZn + (iz-1) * MeshDZ) .lt. Z(e)+errT)) then
            isOnIF = .true.
            end if
         end do

            moviePoints(iP)%center(1) = 0.0
            moviePoints(iP)%center(2) = boundingZn + (iz-1) * MeshDZ
            moviePoints(iP)%normal(1) = 1
            moviePoints(iP)%normal(2) = 1
            moviePoints(iP)%layer = currentLayer
            moviePoints(iP)%isOnInterface = isOnIF
            
           if (verbose >= 2) print*,iP," [", & 
               moviePoints(iP)%center(1),",", & 
               moviePoints(iP)%center(2),"] layer=", moviePoints(iP)%layer
            
            iP = iP + 1
        end do !iz
        
      end subroutine getVideoPoints 
      subroutine getTopography(outpf)
      !Read coordinates of collocation points and fix if there are
      !intersections with inferfaces. Also find normal vectors of segments.
      use GeometryVars
      use gloVars
      use fitting
      use soilVars, only : Z,N,BETA
      use waveNumVars, only : NMAX,vout
      use ploteo10pesos
      use resultVars, only : BouPoints, nBpts, & 
                             mPtfin,bPtini,bPtfin,iPtfin, &
                             ibemMat,trac0vec,IPIVbem 
!     use refSolMatrixVars, only: Bb
      use Gquadrature, only: Gqu_n => Gquad_n, & 
                             Gqu_t => Gqu_t_8, & 
                             Gqu_A => Gqu_A_8
      
      implicit none
      integer, intent(in) :: outpf
      logical :: lexist
      real, dimension(:), allocatable :: auxVector
      real :: nuevoPx,deltaX,deltaZ
      integer :: iXI,iX,e,indice
      real :: XIx, XIy
      character(LEN=100) :: txt
      real     :: errT = 0.001
      logical, dimension(:), allocatable  :: isOnIF
      logical :: huboCambios!, smallEnough
      real, dimension(:), allocatable :: auxA,auxB
      
      real, dimension(2) :: norm_comp
      real, dimension(2) :: ABp !x or y coords of points A and B
      integer :: i,xory,xoryOtro,status
      real :: xA,yA,xB,yB
      real :: interpol  
      
      ! min wavelenght = beta / maxFrec
      huboCambios = .false.
      allocate(auxA(1))
      allocate(auxB(1))
      
      !Read Surface topography
      inquire(file="surface.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="surface.txt",FORM="FORMATTED")
      else
        write(outpf,'(a)')'There is a missing input file. '
        stop 'Check "surface.txt" on Working directory' 
      end if
      
      READ(77,*)
      READ(77,*) nXI !number of points
      ALLOCATE (Xcoord(nXI,2)) !coordinates of points
      ALLOCATE (x(nXI))
      ALLOCATE (y(nXI))
      allocate(auxVector(nXI+1))
      DO iXI = 1,nXI
         READ(77,*) XIx , XIy
         Xcoord(iXI,1) = XIx ; Xcoord(iXI,2) = XIy
      END DO
      close(77)
      
      if (verbose >= 1) then
        write(outpf,'(A,I3,A)') ' We read ',nXI, & 
         ' nodes describing the boundary.'
        DO iXI = 1,nXI
         write(outpf,'(a,F12.8,a,F12.8,a)') & 
         "[",Xcoord(iXI,1),";",Xcoord(iXI,2),"]"
        END DO
      end if
      
      ! Fit Polynomial so we can find more points along the surface.
      x = Xcoord(:,1)
      y = Xcoord(:,2)
      surf_poly = polyfit(x,y,size(x),degree,verbose,outpf)
      !are the polynomial coefficients from lower to top. A0 A1 A2 ... An
      
      ! we need to add a point at every intersection of the 
      ! topography and an interface
      iXI = 1
      DO WHILE (iXI <= nXI-1)!iXI = 1,nXI-1 !para cada segmento
     ! we check if there is a change of medium along segment [iXI,iXI+1]
!        if (verbose >= 2 ) then
!          write(outpf,'(a)')" "
!          write(outpf,'(a)')"next segment "
!        end if
         DO e = 2,size(Z) !en cada estrato
            if (verbose >= 3 ) then
          write(outpf,'(a,F12.8)')"Z =",Z(e)
          write(outpf,'(a,F8.6,a,F8.6,a)')"L: (",y(iXI),",",y(iXI+1),")"
            end if
         
          ! 
          if (anint(y(iXI) * 1000.) == anint(Z(e) * 1000.)) then
           if (verbose >= 2) then; write(outpf,*)"equal"; end if
          else 
      ! si es un segmento que se forma recorriendo los puntos hacia Z+
            if (y(iXI) < Z(e)  .AND. y(iXI+1) > Z(e)) then
               !hay un cambio de estrato
               nuevoPx = splitatY(surf_poly,degree,Z(e),x(iXI),x(iXI+1))
             if (verbose >= 2) then
              write(outpf,'(a,F10.4,a,F10.4,a,F10.4,a,F10.4,a)') & 
              "(dn)Segment [",y(iXI),"-" &
              ,y(iXI+1),"] crosses interface Z= (",nuevoPx,"-",Z(e),")"
             end if
               ! insertamos el nuevo punto en el vector de puntos
               deallocate(auxVector)
               allocate(auxVector(size(x)+1))
               auxVector(1:iXI) = x(1:iXI)
               auxVector(iXI+1) = nuevoPx
               auxVector(iXI+2 : size(auxVector)) = x(iXI+1:size(x))
               deallocate(x)
               allocate(x(size(auxVector)))
               x = auxVector
               
               deallocate(auxVector)
               allocate(auxVector(size(y)+1))
               auxVector(1:iXI) = y(1:iXI)
               auxVector(iXI+1) = polyVal(surf_poly,degree,nuevoPx)
               auxVector(iXI+2 : size(auxVector)) = y(iXI+1:size(y))
               deallocate(y)
               allocate(y(size(auxVector)))
               y = auxVector
               
               nXI = nXI + 1 
            end if 
      ! si es un segmento que se forma recorriendo los puntos hacia Z-
            if (y(iXI) > Z(e)  .AND. y(iXI+1) < Z(e)) then
               !hay un cambio de estrato
               nuevoPx = splitatY(surf_poly,degree,Z(e),x(iXI),x(iXI+1))
             if (verbose >= 2) then
               write(outpf,'(a,F10.4,a,F10.4,a,F10.4,a,F10.4,a)') &
               "(up)Segment [",y(iXI),"-" &
               ,y(iXI+1),"] crosses interface Z= (",nuevoPx,"-",Z(e),")"
             end if
               ! insertamos el nuevo punto en el vector de puntos
               deallocate(auxVector)
               allocate(auxVector(size(x)+1))
               auxVector(1:iXI) = x(1:iXI)
               auxVector(iXI+1) = nuevoPx
               auxVector(iXI+2 : size(auxVector)) = x(iXI+1:size(x))
               deallocate(x)
               allocate(x(size(auxVector)))
               x = auxVector
               
               deallocate(auxVector)
               allocate(auxVector(size(y)+1))
               auxVector(1:iXI) = y(1:iXI)
               auxVector(iXI+1) = polyVal(surf_poly,degree,nuevoPx)
               auxVector(iXI+2 : size(auxVector)) = y(iXI+1:size(y))
               deallocate(y)
               allocate(y(size(auxVector)))
               y = auxVector
               
               nXI = nXI + 1
            end if
          end if
         end do
         iXI = iXI + 1
      END DO
      deallocate(Xcoord)
      nXI = size(x)
      allocate(Xcoord(nXI,2))
      
      ! update the BigSegments nodes
      Xcoord(:,1) = x
      Xcoord(:,2) = y
      
      ! we will not change Xcoord from now on. 
      
      if (verbose >= 1) then
       write(outpf,'(a,I0,a,I0,a)')"We have ",nXI," nodes describing " &
                    ,nXI-1," segments"
       if (Verbose >= 2) then
        DO iXI = 1,nXI
         write(outpf,'(a,F12.8,a,F12.8,a)')"[[",Xcoord(iXI,1),";",Xcoord(iXI,2),"]]"
        END DO
       end if
      end if
      
      deallocate(auxVector)
      
      if (verbose >= 1) then
       CALL chdir("outs",status)
       write(txt,'(a)') 'Geometry.pdf'
       call drawGEOM(txt,.true.,outpf)
       CALL chdir("..")
      end if
      
 ! Unit Normal vector at the midle of BigSegments (shared to subsegments)
      allocate (lengthXI(nXI-1)) !length of Big Segment XI
      allocate (layerXI(nXI-1))
      allocate (midPoint(nXI-1,2)) !center of Big Segment XI
      allocate (isOnIF(nXI-1)) ; isOnIF = .false.
      do iXI = 1,nXI-1
         deltaX = Xcoord(iXI,1) - Xcoord(iXI+1,1)
         deltaZ = Xcoord(iXI,2) - Xcoord(iXI+1,2)
         midPoint(iXI,1) = min(Xcoord(iXI,1),Xcoord(iXI+1,1)) + & 
                           abs(deltaX)/2.
         midPoint(iXI,2) = min(Xcoord(iXI,2),Xcoord(iXI+1,2)) + & 
                           abs(deltaZ)/2.
         lengthXI(iXI) = sqrt(deltaX**2 + deltaZ**2)
         
         e = 0
         do while(Z(e+1) < midPoint(iXI,2))
            e = e + 1
            if (e .gt. N) then
              e = e - 1
              exit
            end if
         end do
         if (e .eq. 0) e = 1
         
         layerXI(iXI) = e
         
         do e=1,N+1
            if((Z(e)-errT .lt. midPoint(iXI,2)) & 
            .and. (midPoint(iXI,2) .lt. Z(e)+errT)) then
            isOnIF(iXI) = .true.
            end if
         end do
      end do
      
      ! rebuild a polynomial that fits Xcoord and the midpoints
      deallocate(x)
      deallocate(y)
      !        nodes + midpoints
      allocate(x(nXI + nXI-1))
      allocate(y(nXI + nXI-1))
      x(1) = Xcoord(1,1)
      y(1) = Xcoord(1,2)
      indice = 1
      do iXI = 2,(2*nXI-1),2
         x(iXI) = midPoint(indice,1)
         x(iXI+1) = Xcoord(indice+1,1)
         y(iXI) = midPoint(indice,2)
         y(iXI+1) = Xcoord(indice+1,2)
         indice = indice + 1
      end do
      surf_poly = polyfit(x,y,size(x),degree,verbose,outpf)
      
!     if (verbose >= 1) then
!      write(outpf,'(a)')'Normal vectors at midpoint of segments [Z+ orientation]'
!     end if
 !normal vector componts at midPoints of BigSegments
      ALLOCATE (normXI(nXI-1,2)) 
      normXI = normalto(midPoint(:,1),midPoint(:,2),nXI-1,surf_poly, & 
               degree,verbose,outpf)
      
      ! Xcoord stores the BigSegment nodes coordinates. 
      ! x and y vectors will be subdivided 
      deallocate(x)
      deallocate(y)
      allocate( x(size(Xcoord(:,1))) )
      allocate( y(size(Xcoord(:,2))) )
      x = Xcoord(:,1)
      y = Xcoord(:,2)
      
      allocate(Nodes_xz_n(size(x),2,2))
      Nodes_xz_n(:,1,1) = x
      Nodes_xz_n(:,2,1) = y
      Nodes_xz_n(1:nXI-1,1,2) = normXI(:,1)
      Nodes_xz_n(1:nXI-1,2,2) = normXI(:,2)
      
      !---
      nBpts = nXI - 1 ! es menos 1 porque se usan los puntos centrales
      bPtini = 1
      if (getInquirePointsSol) bPtini = iPtfin + 1
      if (makeVideo) bPtini = mPtfin + 1
      bPtfin = bPtini + nBpts - 1
      !Boundary points array:
      allocate(BouPoints(nBpts))
      !add center
      BouPoints(:)%center(1) = midPoint(:,1)
      BouPoints(:)%center(2) = midPoint(:,2)
      !borders of segment
      BouPoints(:)%bord_A(1) = Xcoord(1:nXI-1,1)
      BouPoints(:)%bord_A(2) = Xcoord(1:nXI-1,2)
      BouPoints(:)%bord_B(1) = Xcoord(2:nXI,1)
      BouPoints(:)%bord_B(2) = Xcoord(2:nXI,2)
      !add normal
      BouPoints(:)%normal(1) = normXI(:,1)
      BouPoints(:)%normal(2) = normXI(:,2)
      !add length
      BouPoints(:)%length = lengthXI
      !add layer
      BouPoints(:)%layer = layerXI
      BouPoints(:)%isBoundary = .true.
      BouPoints(:)%isOnInterface = isOnIF
      BouPoints(:)%guardarFK = .false.
      BouPoints(:)%guardarMovieSiblings = .false.
      
      do iX=1,nBpts !van vacías porque esto cuenta para cada frecuencia
!       allocate(BouPoints(iX)%FKh(NMAX+1,5)); BouPoints(iX)%FKh = 0
!       allocate(BouPoints(iX)%FKv(NMAX+1,5)); BouPoints(iX)%FKv = 0
        allocate(BouPoints(iX)%FK(1,2*NMAX,5)) 
        allocate(BouPoints(iX)%W(1,5))
        
      ! Coordenadas/pesos de integración Gaussiana:
        allocate(BouPoints(iX)%Gq_xXx_coords(Gqu_n,2))
        allocate(BouPoints(iX)%Gq_xXx_C(Gqu_n))
        
        norm_comp(1)=abs(BouPoints(iX)%bord_B(1)-BouPoints(iX)%bord_A(1)) & 
                      / BouPoints(iX)%length
        norm_comp(2)=abs(BouPoints(iX)%bord_B(2)-BouPoints(iX)%bord_A(2)) & 
                      / BouPoints(iX)%length
        
        if (norm_comp(2) > norm_comp(1)) then
            xory = 2 ! la pendiente es mayormente vertical
            xoryOtro = 1
        else 
            xory = 1 ! la pendiente es mayormente horizontal
            xoryOtro = 2
        end if
        ABp(1) = BouPoints(iX)%bord_A(xory)
        ABp(2) = BouPoints(iX)%bord_B(xory)
        
        do i = 1,Gqu_n !ceros de Legendre (una coordenada):
          BouPoints(iX)%Gq_xXx_coords(i,xory) = (ABp(2)+ABp(1))/2 + &
                                           (ABp(2)-ABp(1))/2 * Gqu_t(i)
                                           
          BouPoints(iX)%Gq_xXx_C(i) = abs(BouPoints(iX)%length)/2 * Gqu_A(i)
        end do
        
        ! la otra coordenada:
        xA = ABp(1)
        yA = BouPoints(iX)%bord_A(xoryOtro)
        xB = ABp(2)
        yB = BouPoints(iX)%bord_B(xoryOtro)
        do i = 1,Gqu_n
      BouPoints(iX)%Gq_xXx_coords(i,xoryOtro) = interpol(xA,yA,xB,yB, &
                                 BouPoints(iX)%Gq_xXx_coords(i,xory))
        end do
        
        if (verbose .ge. 2) then
          write(outpf,'(a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.6,a)') & 
               "[",BouPoints(iX)%bord_A(1),",",BouPoints(iX)%bord_A(2),"]-[", & 
                 BouPoints(iX)%bord_B(1),",",BouPoints(iX)%bord_B(2), & 
                 "] L:", BouPoints(iX)%length, &
                 " n:[",BouPoints(iX)%normal(1),",",BouPoints(iX)%normal(2),"]"
        if (xory .eq. 1) print*,"mayormente horizontal"
        if (xory .eq. 2) print*,"mayormente vertical"          
!       print*,"{",xA,",",yA,"}-{",xB,",",yB,"} Gquad points:"
        do i = 1,Gqu_n
          print*,"Gq",i,"[",BouPoints(iX)%Gq_xXx_coords(i,1), " , ", &
          BouPoints(iX)%Gq_xXx_coords(i,2), "] :: ",BouPoints(iX)%Gq_xXx_C(i)
        end do
        print*,""
        end if !verbose 
        
        ! start saving space for the green functions:
        allocate(BouPoints(iX)%GT_k(nBpts))
        
!       allocate(BouPoints(iX)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); BouPoints(iX)%GT_gq_k = 0
        allocate(BouPoints(iX)%GT_gq(nBpts,Gqu_n,5,2)); BouPoints(iX)%GT_gq = 0
        
      end do !iX
      
      
      deallocate(Vout)
      allocate(Vout(2*nBpts,2))
      allocate(ibemMat(2*nBpts,2*nBpts))
      allocate(trac0vec(2*nBpts))
      allocate(IPIVbem(2*nBpts))
      
      end subroutine getTopography
      subroutine subdivideTopo(iJ,outpf)
      !Read coordinates of collocation points and fix if there are
      !intersections with inferfaces. Also find normal vectors of segments.
      use GeometryVars! x,y,Xcoord,surf_poly,nXI
      use gloVars!, only : verbose, makeVideo, getInquirePointsSol
      use fitting
      use soilVars, only : Z,N,BETA
      use waveNumVars, only : NFREC,DFREC,NMAX,vout
      use resultVars, only : BouPoints, nBpts, & 
                             bPtini,bPtfin,iPtini,iPtfin, &
                             allpoints,nPts,mPtini,mPtfin,NxXxMpts, &
                             ibemMat,trac0vec,IPIVbem
!     use refSolMatrixVars, only: Bb
      use ploteo10pesos
      use Gquadrature, only: Gqu_n => Gquad_n, & 
                             Gqu_t => Gqu_t_8, & 
                             Gqu_A => Gqu_A_8

      implicit none
      integer, intent(in) :: outpf,iJ
      real :: deltaX,deltaZ,maxLen,leniXI
      integer :: J,iXI,iX,e,indice, AllocateStatus,idx
      character(LEN=100) :: txt
      real     :: errT = 0.001
      logical, dimension(:), allocatable  :: isOnIF
      real :: maxFrec
      real, dimension(2) :: midiXI
      logical :: smallEnough, huboCambios
      real, dimension(:), allocatable :: auxA,auxB
      
      real, dimension(2) :: norm_comp
      real, dimension(2) :: ABp !x or y coords of points A and B
      integer :: i,xory,xoryOtro
      real :: xA,yA,xB,yB
      real :: interpol
      integer :: iNstart,iNfinish,iN
      
      huboCambios = .false.
      allocate(isOnIF(1))
      J = iJ
      if (J .lt. 2) J=2
      maxFrec =  DFREC*real(J-1)
      allocate(auxA(1))
      allocate(auxB(1))
      
      !------- subdivide to SEE the smallest wavelenght 
      iXI = 1
      DO WHILE (iXI <= nXI-1) ! for every segment [iXI,iXI+1]
        
           if (verbose >= 3) then
             write(outpf,*)
             write(outpf,'(a,F12.8,a,F12.8,a,F12.8,a,F12.8,a)') & 
                         "at segment: [",x(iXI),",",y(iXI),"] - [", &
                                          x(iXI+1),",",y(iXI+1),"]"
           end if
           !find midPoint of segment
           deltaX = x(iXI) - x(iXI+1)
           deltaZ = y(iXI) - y(iXI+1)
           midiXI(1) = min(x(iXI),x(iXI+1)) + abs(deltaX)/2.
           midiXI(2) = min(y(iXI),y(iXI+1)) + abs(deltaZ)/2.
           e = 0
           do while (Z(e+1) .le. midiXI(2)) 
              e = e + 1
              if (e .gt. N) then
                e = e - 1
                exit
              end if
           end do ! layers
           if (e .eq. 0) e = 1
        ! la longitud máxima para un segmento en este estrato:
           maxLen = BETA(e) / (multSubdiv * maxFrec)
           !safe lock
           if (maxLen < 0.005) then
              write(outpf,'(a,I3)')'e:',e
              write(outpf,'(a,F6.3)')'z:',Z(e)
              write(outpf,'(a,F6.3)')'beta:',BETA(e)
              write(outpf,'(a,F6.3)')'maxlen:',maxLen
              maxLen = max(maxLen,0.005)
              write(outpf,'(a,a)')"  WARNING: Safe lock ative. ", & 
              "maxLen of segments fixed to 0.005m"
           end if
       
         ! aunque los segmentos rectos de integración son del doble de
         ! esta longitud.
        if (verbose >= 2) then
           Write(outpf,'(a,I3,a,I2,a,F10.6,a,F10.6)')"Segmento ",iXI, &
           " está en el estrato: ",e," length=",lengthXI(iXI), &
           " vs max=",maxLen
!       write(outpf,'(a,F10.6,a,F10.6,a)')"[",midiXi(1),",",midiXi(2),"]"
        end if
        
        ! we compare againts segment length
        if (lengthXI(iXI) > maxLen) then
             if (verbose >= 2) then
              write(outpf,'(a,I0,a,F8.2,a,F8.2,a)') & 
              "L(:",iXi,")=",lengthXI(iXI),"> maxLen(",maxLen,"); "
             end if
              !we divide de segment in half and check,
              !if it is not small enough we divide in half again.
              smallEnough = .false.
              indice = 1
              midiXI = 0.0d0
         do while (smallEnough .eqv. .false.)
              !find the mid-point
              deltaX = x(iXI) - x(iXI+1)!; print*,"Dx ",deltaX
              deltaZ = y(iXI) - y(iXI+1)!; print*,"Dz ",deltaZ
              
              midiXI(1) = min(x(iXI),x(iXI+1)) + abs(deltaX)/2.
              midiXI(2) = min(y(iXI),y(iXI+1)) + abs(deltaZ)/2.
              
              leniXI = sqrt(deltaX**2. + deltaZ**2.) / 2.
              
      ! insert new point at the middle of former segment in the list
         deallocate(auxA,auxB, STAT = AllocateStatus)
         allocate(auxA(iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxA = x(1:iXI)
         allocate(auxB(size(x)-iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxB = x(iXI+1:size(x))
         deallocate(x, STAT = AllocateStatus)
         allocate(x(size(auxA)+1+size(auxB)), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         x(1:iXI) = auxA
         x(iXI+1) = midiXI(1)
         x(iXI+2:size(x)) = auxB
         
         deallocate(auxA,auxB, STAT = AllocateStatus)
         allocate(auxA(iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxA = y(1:iXI)
         allocate(auxB(size(y)-iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxB = y(iXI+1:size(y))
         deallocate(y, STAT = AllocateStatus)
         allocate(y(size(auxA)+1+size(auxB)), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         y(1:iXI) = auxA
         y(iXI+1) = midiXI(2)
         y(iXI+2:size(y)) = auxB
         
         deallocate(auxA,auxB, STAT = AllocateStatus)
         allocate(auxA(iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxA = lengthXI(1:iXI)
         allocate(auxB(size(lengthXI)-iXI), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         auxB = lengthXI(iXI+1:size(lengthXI))
         deallocate(lengthXI, STAT = AllocateStatus)
         allocate(lengthXI(size(auxA)+1+size(auxB)), STAT=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "Not enough memory"
         lengthXI(1:iXI) = auxA
         lengthXI(iXI:iXI+1) = leniXI
         lengthXI(iXI+2:size(lengthXI)) = auxB             

              if(verbose >= 2) then
               write(outpf,'(a,F12.8,a,F12.8,a,a,F8.6,a,F8.6)') & 
               "New node at: [",x(iXI+1),",",y(iXI+1),"]" &
                 ," length around it:",lengthXI(iXI),"-", &
                              lengthXI(iXI+1)
              end if
              
         nXI = nXI + 1
         indice = indice + 1
              
              if(verbose >= 2) then
               write(outpf,'(a)')"Now subsegments table looks like:"
               do idx = 1,nXI - 1
               write(outpf,'(a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F8.6)') & 
                 "[",x(idx),",",y(idx),"]-[", & 
                 x(idx+1),",",y(idx+1),"] L:", lengthXI(idx)
               end do
              end if
              !
         if (leniXI <= maxLen) then
            smallEnough = .true.
         else
            if (verbose >= 2) write(outpf,*)"Not small enough"
         end if
        end do !while small enough
             
        huboCambios = .true.
       end if !(lengthXI(iXI) > maxLen)
       
         iXI = iXI + 1 
       end do ! segments
      
      
      if (huboCambios) then
      !-------
      ! encontrar punto centrales de nuevo:
      deallocate(Xcoord)
      deallocate(lengthXI)
      deallocate(layerXI)
      deallocate(midPoint)
      deallocate(isOnIF)
      
      nXI = size(x)
      allocate(Xcoord(nXI,2))
      
      ! update the BigSegments nodes
      Xcoord(:,1) = x
      Xcoord(:,2) = y
      
      allocate (lengthXI(nXI-1)) !length of Big Segment XI
      allocate (layerXI(nXI-1))
      allocate (midPoint(nXI-1,2)) !center of Big Segment XI
      allocate (isOnIF(nXI-1)) ; isOnIF = .false.
      !     print*,"Length is found."
      do iXI = 1,nXI-1
         deltaX = Xcoord(iXI,1) - Xcoord(iXI+1,1)
         deltaZ = Xcoord(iXI,2) - Xcoord(iXI+1,2)
         midPoint(iXI,1) = min(Xcoord(iXI,1),Xcoord(iXI+1,1)) + & 
                           abs(deltaX)/2.
         midPoint(iXI,2) = min(Xcoord(iXI,2),Xcoord(iXI+1,2)) + & 
                           abs(deltaZ)/2.
         lengthXI(iXI) = sqrt(deltaX**2 + deltaZ**2)
         
         e = 0
         do while(Z(e+1) < midPoint(iXI,2))
            e = e + 1
            if (e .gt. N) then
               e = e - 1
               exit
            end if
         end do
         if (e .eq. 0) e = 1
         layerXI(iXI) = e
         
         do e=1,N+1
            if((Z(e)-errT .lt. midPoint(iXI,2)) & 
            .and. (midPoint(iXI,2) .lt. Z(e)+errT)) then
            isOnIF(iXI) = .true.
            end if
         end do
      end do
      
      ! rebuild a polynomial that fits Xcoord and the midpoints
      deallocate(x)
      deallocate(y)
      !        nodes + midpoints
      allocate(x(nXI + nXI-1))
      allocate(y(nXI + nXI-1))
      x(1) = Xcoord(1,1)
      y(1) = Xcoord(1,2)
      indice = 1
      do iXI = 2,(2*nXI-1),2
         x(iXI) = midPoint(indice,1)
         x(iXI+1) = Xcoord(indice+1,1)
         y(iXI) = midPoint(indice,2)
         y(iXI+1) = Xcoord(indice+1,2)
         indice = indice + 1
      end do
      surf_poly = polyfit(x,y,size(x),degree,verbose,outpf)

      !normal vector componts at midPoints of BigSegments
      deallocate(normXI)
      ALLOCATE (normXI(nXI-1,2)) 
!     normXI = normalto(midPoint(:,1),midPoint(:,2),nXI-1,surf_poly, & 
!              degree,verbose,outpf)
      
      ! en lugar de hacer esto, vamos a heredar las normales para cada rango de ordenadas --------------------
      ! Nodes_xz_n
      iNstart = 1
      iNfinish = size(Nodes_xz_n(:,1,1)) - 1
      do iX = 1,size(midPoint(:,1))
      do iN = iNstart,iNfinish
      if(Nodes_xz_n(iN,1,1) .le. midPoint(ix,1) .and. midPoint(ix,1) .le. Nodes_xz_n(iN+1,1,1)) then
      if (verbose .ge. 2) print*,midPoint(ix,1)," entre {",Nodes_xz_n(iN,1,1)," - ",Nodes_xz_n(iN+1,1,1),"}"
      normXI(iX,1) = Nodes_xz_n(iN,1,2)
      normXI(iX,2) = Nodes_xz_n(iN,2,2)
      end if
      end do
      end do
      
      if (verbose >= 1) then
       CALL chdir("outs")
       write(txt,'(a,I0,a)') 'Division_at[',iJ,'Hz].pdf'
       call drawGEOM(txt,.false.,outpf)
       CALL chdir("..")
        write(outpf,'(a,I0,a)')"Boundary was segmented, now it is made of: (",nXI-1,") elements"
      end if!
      if (verbose >= 2) then
        do idx = 1,nXI - 1
          write(outpf,'(a,F7.3,a,F7.3,a,F7.3,a,F7.3,a,F6.2,a,F7.3,a,F7.3,a)') & 
            "[", Xcoord(idx,1),",", Xcoord(idx,2),"]->[", & 
               Xcoord(idx+1,1),",", Xcoord(idx+1,2),"] L:", &
               lengthXI(idx)," mid: [",midPoint(idx,1),",",midPoint(idx,2),"]"
        end do
      end if
      deallocate(x)
      deallocate(y)
      !        nodes + midpoints
      allocate(x(nXI))
      allocate(y(nXI))
      x = Xcoord(:,1)
      y = Xcoord(:,2)    
      !---
      nBpts = nXI - 1 ! es menos 1 porque se usan los puntos centrales
      bPtini = 1
      if (getInquirePointsSol) bPtini = iPtfin + 1
      if (makeVideo) bPtini = mPtfin + 1
      bPtfin = bPtini + nBpts - 1
      !Boundary points array:
      do iX = 1,size(BouPoints) ! para no chorrear memoria
!       deallocate(BouPoints(iX)%FKh)
!       deallocate(BouPoints(iX)%FKv)
        deallocate(BouPoints(iX)%FK)
        deallocate(BouPoints(iX)%W)
        deallocate(BouPoints(iX)%Gq_xXx_coords)
        deallocate(BouPoints(iX)%Gq_xXx_C)
!       deallocate(BouPoints(iX)%GT_gq_k)
        deallocate(BouPoints(iX)%GT_gq)
!       do iXi=1,size(BouPoints(iX)%GT_k)
!         if(allocated(BouPoints(iX)%GT_k(ixi)%qlmk)) deallocate(BouPoints(iX)%GT_k(iXi)%qlmk)
!       end do
        deallocate(BouPoints(iX)%GT_k)
      end do
      deallocate(BouPoints)
      allocate(BouPoints(nBpts))

      !add center
      BouPoints(:)%center(1) = midPoint(:,1)
      BouPoints(:)%center(2) = midPoint(:,2)
      !borders of segment
      BouPoints(:)%bord_A(1) = Xcoord(1:nXI-1,1)
      BouPoints(:)%bord_A(2) = Xcoord(1:nXI-1,2)
      BouPoints(:)%bord_B(1) = Xcoord(2:nXI,1)
      BouPoints(:)%bord_B(2) = Xcoord(2:nXI,2)
      !add normal
      BouPoints(:)%normal(1) = normXI(:,1)
      BouPoints(:)%normal(2) = normXI(:,2)
      !add length
      BouPoints(:)%length = lengthXI
      !add layer
      BouPoints(:)%layer = layerXI
      BouPoints(:)%isBoundary = .true.
      BouPoints(:)%isOnInterface = isOnIF
      BouPoints(:)%guardarFK = .false.
      BouPoints(:)%guardarMovieSiblings = .false.
      
      ! para resolver las densidades de fuerza IBEM:
      do iX=1,nBpts !van vacías porque esto cuenta para cada frecuencia
        ! start saving space for the green functions:
        allocate(BouPoints(iX)%GT_k(nBpts))
!       allocate(BouPoints(iX)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); BouPoints(iX)%GT_gq_k = 0
        allocate(BouPoints(iX)%GT_gq(nBpts,Gqu_n,5,2)); BouPoints(iX)%GT_gq = 0
          
!       allocate(BouPoints(iX)%FKh(NMAX+1,5)); BouPoints(iX)%FKh = 0
!       allocate(BouPoints(iX)%FKv(NMAX+1,5)); BouPoints(iX)%FKv = 0
        allocate(BouPoints(iX)%FK(1,2*NMAX,5))  
        allocate(BouPoints(iX)%W(1,5))
        
        allocate(BouPoints(iX)%Gq_xXx_coords(Gqu_n,2))
        allocate(BouPoints(iX)%Gq_xXx_C(Gqu_n))
       
      ! Coordenadas de los puntos de integración Gaussiana.
        norm_comp(1)=abs(BouPoints(iX)%bord_B(1)-BouPoints(iX)%bord_A(1)) & 
                      / BouPoints(iX)%length
        norm_comp(2)=abs(BouPoints(iX)%bord_B(2)-BouPoints(iX)%bord_A(2)) & 
                      / BouPoints(iX)%length
        
        if (norm_comp(2) > norm_comp(1)) then
  !           print*," la pendiente es mayormente vertical "
            xory = 2 ! la pendiente es mayormente vertical
            xoryOtro = 1
        else 
  !           print*," la pendiente es mayormente horizontal "
            xory = 1 ! la pendiente es mayormente horizontal
            xoryOtro = 2
        end if
        ABp(1) = BouPoints(iX)%bord_A(xory)
        ABp(2) = BouPoints(iX)%bord_B(xory)
        
        do i = 1,Gqu_n !ceros de Legendre (una coordenada):
          BouPoints(iX)%Gq_xXx_coords(i,xory) = (ABp(2)+ABp(1))/2 + &
                                           (ABp(2)-ABp(1))/2 * Gqu_t(i)
                                           
          BouPoints(iX)%Gq_xXx_C(i) = abs(BouPoints(iX)%length)/2 * Gqu_A(i)
        end do
        
        ! la otra coordenada:
        xA = ABp(1)
        yA = BouPoints(iX)%bord_A(xoryOtro)
        xB = ABp(2)
        yB = BouPoints(iX)%bord_B(xoryOtro)
        do i = 1,Gqu_n
      BouPoints(iX)%Gq_xXx_coords(i,xoryOtro) = interpol(xA,yA,xB,yB, &
                                 BouPoints(iX)%Gq_xXx_coords(i,xory))
        end do
        
        if (verbose .ge. 2) then
        write(outpf,'(a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.8,a,F12.6,a)') & 
               "[",BouPoints(iX)%bord_A(1),",",BouPoints(iX)%bord_A(2),"]-[", & 
                 BouPoints(iX)%bord_B(1),",",BouPoints(iX)%bord_B(2), & 
                 "] L:", BouPoints(iX)%length, &
                 " n:[",BouPoints(iX)%normal(1),",",BouPoints(iX)%normal(2),"]"
        if (xory .eq. 1) print*,"mayormente horizontal"
        if (xory .eq. 2) print*,"mayormente vertical"          
!       print*,"{",xA,",",yA,"}-{",xB,",",yB,"} Gquad points:"
        do i = 1,Gqu_n
          print*,"Gq",i,"[",BouPoints(iX)%Gq_xXx_coords(i,1), " , ", &
          BouPoints(iX)%Gq_xXx_coords(i,2), "] :: ",BouPoints(iX)%Gq_xXx_C(i)
        end do
        print*,""
        end if
      end do !iX
      
      ! G para resolver el campo difractado por topografía
      do iX=1,nPts
        if(allocated(allpoints(ix)%GT_k)) deallocate(allpoints(ix)%GT_k) 
        allocate(allpoints(ix)%GT_k(nBpts))
      
!       deallocate(allpoints(ix)%GT_gq_k)
!         allocate(allpoints(ix)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); allpoints(ix)%GT_gq_k = 0
      end do!
      do iX=iPtini,iPtfin 
        deallocate(allpoints(ix)%GT_gq)
          allocate(allpoints(ix)%GT_gq(nBpts,Gqu_n,5,2)); allpoints(ix)%GT_gq = 0
      end do
      
      if (makeVideo) then
        do iX = mPtini,mPtfin
          deallocate(allpoints(ix)%GT_gq_mov)
          allocate(allpoints(ix)%GT_gq_mov(nBpts,Gqu_n,5,2,NxXxMpts)); allpoints(ix)%GT_gq_mov = 0
        end do
      end if
      
      !también en necesario actulizar el tamaño del vector
      !de términos independientes
      
      deallocate(ibemMat, trac0vec, IPIVbem,Vout)
      allocate(Vout(2*nbpts,2))
      allocate(ibemMat(2*nBpts,2*nBpts))
      allocate(trac0vec(2*nBpts))
      allocate(IPIVbem(2*nBpts))
      end if ! huboCambios
      
      end subroutine subdivideTopo
      subroutine allocintegPoints(iJ)
      use soilVars, only : BETA
      use resultVars, only : allpoints,BouPoints,npts,nBpts
      use waveNumVars, only : NMAX,DFREC
      use Gquadrature, only: WLmulti,Gqu_n => Gquad_n
      
      implicit none
      integer, intent(in) :: iJ
      integer :: J,i_X,iXi
      real :: MinWaveLenght,distance
      J = iJ
      if (J .lt. 2) J = 2
      
!     if(.not. allocated(BouPoints(1)%GT_k(1)%qlmk)) then
      do iXi = 1,nBpts
      do i_x = 1,nBpts
      !Si la distancia entre segmentos es menor a CERCA usamos la quadratura
         
        MinWaveLenght = min(BETA(BouPoints(i_X)%layer),BETA(BouPoints(iXi)%layer)) / (DFREC*real(J-1))
        distance = sqrt((BouPoints(i_X)%center(1) - BouPoints(iXi)%center(1))**2 + &
                        (BouPoints(i_X)%center(2) - BouPoints(iXi)%center(2))**2)
        
        if (distance .lt. MinWaveLenght * WLmulti) then !use quadratura
!         if(.not. allocated(BouPoints(i_X)%GT_k(iXi)%qlmk)) &
!         allocate(BouPoints(i_X)%GT_k(iXi)%qlmk(Gqu_n,5,2,nmax+1))
          BouPoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .true.
!         allocate(BouPoints(i_X)%GT_gq(nBpts,Gqu_n,5,2))
        else !use integración directa
!         if(.not. allocated(BouPoints(i_X)%GT_k(iXi)%qlmk)) &
!         allocate(BouPoints(i_X)%GT_k(iXi)%qlmk(1,5,2,nmax+1))
          BouPoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .false.
        end if
      end do !i_x
      
      do i_x = 1,npts
      !Si la distancia entre segmentos es menor a CERCA usamos la quadratura
         
        MinWaveLenght = min(BETA(allpoints(i_X)%layer),BETA(BouPoints(iXi)%layer)) / (DFREC*real(J-1))
        distance = sqrt((allpoints(i_X)%center(1) - BouPoints(iXi)%center(1))**2 + &
                        (allpoints(i_X)%center(2) - BouPoints(iXi)%center(2))**2)
        
        if (distance .lt. MinWaveLenght * WLmulti) then !use quadratura
!         if(.not. allocated(allpoints(i_X)%GT_k(iXi)%qlmk)) &
!         allocate(allpoints(i_X)%GT_k(iXi)%qlmk(Gqu_n,5,2,nmax+1))
          allpoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .true.
!         allocate(allpoints(i_X)%GT_gq(nBpts,Gqu_n,5,2))
        else !use integración directa
!         if(.not. allocated(allpoints(i_X)%GT_k(iXi)%qlmk)) &
!         allocate(allpoints(i_X)%GT_k(iXi)%qlmk(1,5,2,nmax+1))
          allpoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .false.
        end if
      end do !i_x
      
      end do !ixi
!     end if
      end subroutine allocintegPoints
      function interpol(xA,yA,xB,yB,x)
      implicit none
      real :: interpol
      real, intent(in) :: xA,yA,xB,yB,x
      real :: m
      !interpolación lineal de la ordenada de x que está entre A y B
      m = (yB-yA) / (xB-xA)
      interpol = yA + m * (x-xA)
      end  function interpol
      subroutine waveAmplitude(outpf)                                   !    FIX NFREC
      use wavelets !las funciones: ricker, fork
      use waveVars, only : Dt,Uo
      use waveNumVars, only : NFREC,DFREC
      use gloVars, only : verbose
      use ploteo10pesos
      implicit none
      integer, intent(in) :: outpf
      integer :: i
      character(LEN=9)                             :: xAx,yAx,logflag
      character(LEN=100)                           :: titleN
!     real :: factor
      ! Amplitude of incident wave
      ! prepare the signal we will use as incident wave amplitude.
      ! el tamaño del ricker es 2*NFREC porque incluirá el conjugado
      if (verbose >= 1) then
       write(outpf,'(a)')' '
       write(outpf,'(a)')'Incident wave amplitude: '
      end if
      
      call ricker(NFREC*2,outpf) ! the ricker wavelet saved on Uo
      
      if (verbose >= 1) then
      CALL chdir("outs")
!      OPEN(31,FILE="rick_time.txt",FORM="FORMATTED")
!      write(31,'(I2)') 1
!      write(31,'(a)') "amplitude"
!      write(31,'(F15.8)') Dt
!      do i = 1,size(Uo)
!         write(31,'(F50.16,2x,F50.16)') real(Uo(i)),aimag(Uo(i))
!      end do       
!      close (31)
       ! plot it to a pdf file:
       write(titleN,'(a)') 'rick_time.pdf'
       xAx = 'time[sec]'
       yAx = 'amplitude'
       call plotXYcomp(Uo,Dt,size(Uo),titleN,xAx,yAx,1200,800)
       CALL chdir("..")
!      CALL SYSTEM ('../plotXYcomp rick_time.txt rick_time.pdf time[sec] amplitude 1200 800')
      end if
      
      call fork(size(Uo),Uo,-1,verbose,outpf) ! fft into frequency 
      
      if (verbose >= 1) then
       CALL chdir("outs")
       OPEN(32,FILE="rick_frec.txt",FORM="FORMATTED")
       write(32,'(I2)') 1
       write(32,'(a)') "amplitude"
       write(32,'(F15.8)') DFREC  !1/(Dt * size(Uo))
       do i = 1, size(Uo)
          write(32,'(F50.16,2x,F50.16)') real(Uo(i)),aimag(Uo(i))
       end do
       close (32)
       ! plot it too
       write(titleN,'(a)') 'rick_frec.pdf'
       xAx = 'frec[Hz] '
       yAx = 'amplitude'
       logflag = 'logx     '
       call plotSpectrum(Uo,DFREC,size(Uo),int(size(Uo)/2),titleN,xAx,yAx,logflag,1200,800)
       CALL chdir("..")
!      logflag = 'logy     '
!      call plotSpectrum(Uo,DFREC,size(Uo),int(size(Uo)/2),titleN,xAx,yAx,logflag,1200,800)
                        
!      call SYSTEM('../plotSpectrum rick_frec.txt rick_frec.pdf frecuency[Hz] amplitude logx 1200 800')
!      !CALL SYSTEM ('../plotXYcomp rick_frec.txt & 
!      !            rick_frec_ugly.pdf time[sec] amplitude 1200 800')
      end if
      
!     ! test............................................................
!     ! before anything, we have to beautify the signal
!     n = size(Uo)
!     factor = sqrt(real(n))
!     Uo(1) = Uo(1)*factor
!     do i=2,n/2+1
!       Uo(i) = Uo(i)*factor! * Gij(i-1)
!       Uo(n-i+2) = conjg(Uo(i))
!     end do!
!     do i = 1, size(Uo)
!         write(6,'(F50.16,2x,F50.16)') real(Uo(i)),aimag(Uo(i))
!      end do
!     if (verbose >= 3) then
!     call fork(size(Uo), Uo,+1,verbose,outpf) !!fft back
!     Uo(:) = Uo(:)/factor
!     ! aquí se quitaría el efecto de la frecuencia imaginaria...
!     
!      CALL chdir("outs")
!      OPEN(33,FILE="rick_time_back.txt",FORM="FORMATTED")
!      write(33,'(I2)') 1
!      write(33,'(a)') "amplitude"
!      write(33,'(F15.8)') Dt
!      do i = 1,size(Uo)
!       write(33,'(F50.16,2x,F50.16)') real(Uo(i)),aimag(Uo(i))
!      end do
!      close (33)
!      ! plot it to an pdf file:
!      write(titleN,'(a)') 'rick_timB.pdf'
!      xAx = 'time[sec]'
!      yAx = 'amplitude'
!      call plotXYcomp(Uo,Dt,size(Uo),titleN,xAx,yAx,1200,800)
!     end if
!     stop 0
!     ! ................................................................
      
      if (verbose >= 1) then
       write(outpf,'(a)') ""
       write(outpf,'(A,I4)') ' Ricker(w) signal size is :',size(Uo)
      end if
      end subroutine waveAmplitude
      subroutine stratumDiffField (iPxi, &
                                   PX, nPX, task, &
                                   J, cOME, outpf, V)
      ! task = -1
      ! Toma una fuente XI y devuelve un vector V_[2*nBpts,2] de tracciones sobre puntos X
      ! 
      ! task =  0
      ! Almacena el campo diffractado (estratigrafía) en W y WmovieSiblings
      
      use gloVars, only: verbose,PI
      use resultvars, only:nBpts,MecaElem,Punto,Boupoints,Hf,Hk, NxXxMpts
      use sourceVars
      use Gquadrature, only: Gquad_n
      use waveNumVars, only : NMAX,DK,expK
      use refSolMatrixVars, only : B,IPIV,Ak
      use soilVars, only: N
      use meshVars, only: MeshDXmultiplo,nx
      use wavelets !fork
      implicit none
      integer, intent(in) :: iPxi,nPX,task,J,outpf
      complex*16, intent(in)  :: cOME
      type(Punto), dimension(nPX), intent(inout) :: PX
      integer, dimension(5,2) :: sign
      integer :: renglon,iP_x,ef,GQmx,iGq,dir,ik,iMec,xXx,i
      logical :: intf
      complex*16, dimension(2*nBpts,2) :: V
      complex*16, dimension(2*nmax,5) :: auxK
      real :: zf,xf,k,factor,peso
      real, dimension(2) :: nf
      type(MecaElem)  :: calcMecaAt_xi_zi, Meca_diff, calcFF, Meca_FF
      !  fza horiz     fza vert
      sign(1,1) = -1 ; sign(1,2) = 1  !W        impar           par
      sign(2,1) = 1  ; sign(2,2) = -1 !U          par           impar
      sign(3,1) = -1 ; sign(3,2) = 1  !s33      impar           par
      sign(4,1) = 1  ; sign(4,2) = -1 !s31        par           impar
      sign(5,1) = -1 ; sign(5,2) = 1  !s11      impar           par
      factor = sqrt(2.*NMAX*2)
      renglon = 1
      V = 0 
      do iP_x = 1,nPX ! para cada punto X
      
      if (iPxi .eq. iP_x) then
        V(renglon+0,1) = 0.5; V(renglon+0,2) = 0.0
        V(renglon+1,1) = 0.0; V(renglon+1,2) = 0.5
        renglon = renglon + 2
        cycle 
      end if ! en caso de que fueran iguales, la solución es analítica
      
      if (iPxi .eq. 0) then
        ! se trata de la fuente real
        zf = zfsource; xf = xfsource; ef = efsource; intf = intfsource
        nf(2) = nzfsource; nf(1) = nxfsource
        GQmx = 1
      else
        ! se trata de una fuente virtual con índice iPxi
         ! ¿Usar quadratura Gaussiana entre PX y PXI?
        if (PX(iP_x)%GT_k(iPxi)%shouldUseQuadrature) then
          zf = -1 !indicamos que la fuente es cada punto de la quadratura
          xf = -1
        GQmx = Gquad_n
        else
          zf = Boupoints(iPxi)%center(2)
          xf = Boupoints(iPxi)%center(1)
        GQmx = 1
        end if
          ef = BouPoints(iPxi)%layer
        intf = Boupoints(iPxi)%isOnInterface
         nf(2) = Boupoints(iPxi)%normal(2)
         nf(1) = Boupoints(iPxi)%normal(1)
      end if
      
      !campo difractado por la estratigrafía
      do dir = 1,2 ! componete de Fuerza horizontal; vertical
      do iGq = 1,GQmx
        if (GQmx .gt. 1) then !cada punto de la quadratura
          zf = BouPoints(iPxi)%Gq_xXx_coords(iGq,2)
          xf = BouPoints(iPxi)%Gq_xXx_coords(iGq,1)
        end if
     
        auxK = 0
      Do ik=1,nmax+1 !WAVE NUMBER LOOP-
         k = real(ik-1,8) * dK
         if (ik .eq. 1) k = dk * 0.001
        
        call vectorB_force(B,zf,ef,intf,dir,cOME,k)  
        B = matmul(AK(:,:,ik),B)
        Meca_diff = calcMecaAt_xi_zi(B,&  
                     PX(iP_x)%center(2),& 
                     PX(iP_x)%layer,cOME,k,outpf)
        ! campo libre 
        Meca_ff = calcFF(zf,ef,dir, & 
                     PX(iP_x)%center(2),& 
                     PX(iP_x)%layer,k,cOME)
        
        auxK(ik,1:5) =  Meca_ff%Rw(1:5) + Meca_diff%Rw(1:5)
      end do !ik
         ! completar la crepa
      !        calculado              crepa 
      !   _________|____________ _______|________  
      !    1 2 3 ... NMAX NMAX+1 0 0 0 0 0 2*NMAX     
      ! k= 0 1 2 ... N/2-1 -N/2 ... -3 -2 -1
      do iMec = 1,5
         auxK(nmax+1:nmax*2,iMec) = sign(iMec,dir)* auxK(nmax+1:2:-1,iMec)
         auxK(:,iMec) = auxK(:,iMec) * expK ** xf / expK ** PX(iP_x)%center(1) * cmplx(Hf(J) * Hk,0.0,8)
         !                             info X-coord        info X-coord           taper
      end do !iMec
        
        ! guardar el FK si es que es interesante verlo
      if (PX(iP_x)%guardarFK) then
        PX(iP_x)%FK(J,1:nmax,1:5) = auxK(1:nmax,1:5)* nf(dir) + & 
        PX(iP_x)%FK(J,1:nmax,1:5)
      end if
        
        ! K -> X
        auxK = auxK * factor
      do iMec = 1,5
        call fork(2*nmax,auxK(:,iMec),+1,verbose,outpf)
      end do
        auxK = auxK / factor
        
      if (task .eq. -1) then ! obtener vector de tracciones
       !  m=0 m=1 :: dirección de aplicación de la fuerza :: dir
       ! V( 1  2 ) l=0 |__  
       !  ( 1  2 ) l=1 |   X=1
       !  ( 3  4 ) l=0  |__ 
       !  ( 3  4 ) l=1  |   X=2
       !  ( .  . ) ... componente de la tracción en X
       !  ( .  . )
       !   ...
       
       if (GQmx .eq. 1) then
         peso = BouPoints(iPxi)%length
       else
         peso = BouPoints(iPxi)%Gq_xXx_C(iGq)
       end if
       ! componente horizontal de la tracción en X
       V(renglon+0,dir) = V(renglon+0,dir) + & 
        (auxK(1,5-0) * BouPoints(iP_x)%normal(1) + &
         auxK(1,4-0) * BouPoints(iP_x)%normal(2)) * peso
       ! componente vertical
       V(renglon+1,dir) = V(renglon+1,dir) + & 
        (auxK(1,5-1) * BouPoints(iP_x)%normal(1) + &
         auxK(1,4-1) * BouPoints(iP_x)%normal(2)) * peso
         ! componente horizontal                   vertical
         ! t(n)  =  nx Sxx + nz Szx    |  t(n)  =  nx Sxz + nz Szz
         !    i=x                      |     i=z
         !---------------------------------------------------------
         !    i -> x :: Sxx   Szx      |     i -> z :: Sxz   Szz
         !               5     4       |                4     3
         !      0       5-0   4-0      |         1     5-1   4-1
         !---------------------------------------------------------
         
!     !   ,--- componente de la tracción : l
!     !   |,--- (1),(3) dirección de la fuerza : m
!     !   ||    ,--- cara
!     !   ||    |,--- fza
!     !  Txx = Sxx1 nx + Szx1 nz  |___ (1) fza horizontal  .
!     !  Tzx = Szx1 nx + Szz1 nz  |                        . 
!     !  Txz = Sxx3 nx + Szx3 nz |___ (3) fza vertical     .
!     !  Tzz = Szx3 nx + Szz3 nz |                         . 
!     !
!     !  T_lm = S_lkm * n_k
!     !       = S_l1m * n_1 + S_l3m * n3
!     !         __|__         __|__
!     !      s11     s31   s13    s33
!     !      s11     s31   s31    s33  (son equivalentes)
!     !       5       4     4      3   (indice en GT_gq(_,_,i,_) )
!     !       0       1     0      1   (indice de submatriz: l )
      else !task = 0
        ! guardamos el campo difractado
        if(PX(iP_x)%guardarMovieSiblings) then
          ! (1) reordenar la crepa existente en X
          auxK = cshift(auxK,SHIFT=nmax+1,DIM=1)
          ! (2) guardar sólo los interesantes
          xXx = 1
          do i=nmax-(nxXxmpts-1)/2* meshdxmultiplo,nmax+(nxXxmpts-1)/2* meshdxmultiplo,meshdxmultiplo
            if (iPxi .eq. 0) then !fuente real
              if (dir .eq. 1) PX(iP_x)%WmovieSiblings(xXx,J,:) = 0
              PX(iP_x)%WmovieSiblings(xXx,J,1:5) = auxK(i,1:5) * nf(dir) + &
              PX(iP_x)%WmovieSiblings(xXx,J,1:5)
            else ! func de green (sólo desplazamientos)
              PX(iP_x)%GT_gq_mov(iPxi,iGq,1:5,dir,xXx) = auxK(i,1:5) * nf(dir) !+ &
!             PX(iP_x)%GT_gq_mov(iPxi,iGq,1:2,dir,xXx) 
            end if
            xXx = xXx + 1
          end do 
        else !es un inquirePoint
          if (iPxi .eq. 0) then !fuente real
            if (dir .eq. 1) PX(iP_x)%W(J,:) = 0
            PX(iP_x)%W(J,1:5) = auxK(1,1:5) * nf(dir) + &
            PX(iP_x)%W(J,1:5)
          else ! func de green
            PX(iP_x)%GT_gq(iPxi,iGq,1:5,dir) = auxK(1,1:5) * nf(dir) !+ &
!           PX(iP_x)%GT_gq(iPxi,iGq,1:2,dir) 
          end if
        end if
      end if ! task
      end do !dir
      end do !iGq
       renglon = renglon + 2
      end do !iP_x
      end subroutine stratumDiffField
!     subroutine FreeField(z_f,e_f,iP_xi,iGq,P,Npuntos_X,cOME)
!     !                               |
!     use soilVars !                  `--- (0): Fuente real
!     use gloVars  !                   (iP_xi): Fuente virtual
!     use waveNumVars, only : NFREC,NMAX,DK
!     use resultvars, only: Punto
!     use Gquadrature, only: Gquad_n
!     
!     implicit none
!     integer,    intent(in)     :: e_f,iP_xi,iGq,Npuntos_X
!     real,       intent(in)     :: z_f
!     complex*16, intent(in)     :: cOME
!     type(Punto), dimension(Npuntos_X), intent(inout)  :: P
!     
!     integer                    :: ik,iP
!     real                       :: k
!     real, dimension(Npuntos_X) :: z_i,SGNz
!     complex*16, dimension(NMAX+1,Npuntos_X) :: G11,G33,G31, & 
!                                         S333,S313,S113,s331,s131,s111
!     complex*16, dimension(Npuntos_X) :: egamz,enuz
!     real                       :: L2M
!     complex*16                 :: gamma,nu,DEN
!     complex*16                 :: omeAlf,omeBet
!     
!     G11=0;G31=0;G33=0
!     S333=0;S313=0;S113=0;s331=0;s131=0;s111=0
!     egamz=0;enuz=0;z_i=0;SGNz=0
!     egamz=cmplx(1.0,0.0,8);enuz=cmplx(1.0,0.0,8)
!     
!     DEN = 4.0*PI*RHO(e_f)*cOME**2.0
!     omeAlf = cOME**2/ALFA(e_f)**2.0
!     omeBet = cOME**2/BETA(e_f)**2.0
!     L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
!     
!       where (P(:)%layer .eq. e_f)
!          z_i = P(:)%center(2) - z_f
!       end where
!     !
!     where (z_i .ne. 0)
!       SGNz = z_i / ABS(z_i)
!     end where!
!     
!     !loop en k para llenar el plano FK
!     do ik =1,nmax+1
!        k = real(ik-1,8) * dK
!        if (ik .eq. 1) k = dk * 0.001
!     
!       gamma = sqrt(omeAlf - k**2.0)
!       nu = sqrt(omeBet - k**2.0)
!       if(aimag(gamma).gt.0.0)gamma = -gamma
!       if(aimag(nu).gt.0.0)nu=-nu
!               
!       where (P(:)%layer .eq. e_f) 
!         egamz = exp(-UI*gamma*ABS(z_i))
!         enuz = exp(-UI*nu*ABS(z_i)) 
!          
!         G11(ik,:) = -UI/DEN * (k**2.0/gamma*egamz + nu*enuz) 
!         G33(ik,:) = -UI/DEN * (gamma*egamz + k**2.0/nu*enuz) 
!         G31(ik,:) = -UI/DEN * SGNz*k*(egamz - enuz)
!        S333(ik,:) = -UR/DEN * SGNz*( &
!                   (gamma**2.0*L2M + k**2.0*lambda(e_f))* egamz &
!                 + (2.0*amu(e_f)*k**2.0)* enuz &
!                   ) 
!        S313(ik,:) = -UR/DEN * amu(e_f) * ( &
!                   (2.0*k*gamma)* egamz &
!                 - (k/nu*(nu**2.0-k**2.0))* enuz &
!                   ) 
!        S113(ik,:) = -UR/DEN * SGNz * ( &
!                   (k**2.0*L2M + gamma**2.0*lambda(e_f))* egamz &
!                 + (-2.0*amu(e_f)*k**2.0)* enuz &
!                   )
!        s331(ik,:) = -UR/DEN * ( &
!                   (k*gamma*L2M + lambda(e_f)*k**3.0/gamma)* egamz &
!                 + (-2.0*amu(e_f)*k*nu)* enuz &
!                   )
!        s131(ik,:) = -UR/DEN * amu(e_f)*SGNz * ( &             
!                   (2.0*k**2.0)* egamz &
!                   + (nu**2.0-k**2.0)* enuz &
!                   )
!        s111(ik,:) = -UR/DEN * ( & 
!                   (k*gamma*lambda(e_f)+L2M*k**3.0/gamma)* egamz &
!                   + (2.0*amu(e_f)*k*nu)* enuz & 
!                   )
!       end where
!       
!     end do !ik
!    
!     if (iP_xi .eq. 0) then !the case when we the source is the source 
!      do iP = 1,Npuntos_X !Npts
!       P(iP)%FKv(:,1) = G33(:,iP)
!       P(iP)%FKv(:,2) = G31(:,iP)
!       P(iP)%FKv(:,3) = S333(:,iP)
!       P(iP)%FKv(:,4) = S313(:,iP)
!       P(iP)%FKv(:,5) = S113(:,iP)
!       
!       P(iP)%FKh(:,1) = G31(:,iP)
!       P(iP)%FKh(:,2) = G11(:,iP)
!       P(iP)%FKh(:,3) = S331(:,iP)
!       P(iP)%FKh(:,4) = S131(:,iP)
!       P(iP)%FKh(:,5) = S111(:,iP)
!      end do
!     else ! the case when the source is at some point iP_xi 
!      do iP = 1,Npuntos_X
!         !vertical
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,1,2,:) = G33(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,2,2,:) = G31(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,3,2,:) = S333(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,4,2,:) = S313(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,5,2,:) = S113(:,iP)
!         !horizontal
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,1,1,:) = G31(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,2,1,:) = G11(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,3,1,:) = S331(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,4,1,:) = S131(:,iP)
!          P(iP)%GT_k(iP_xi)%qlmk(iGq,5,1,:) = S111(:,iP)
!        
!      end do
!     end if!
!     
!     end subroutine FreeField
      
      ! P plana incidente
!     print*,"Zstart= ",Zstart
!     stop 0
      ! for every element of the geometry subdivided at this frec. 
      ! and other iteresting points
!     do i = 1,allpoints(nJ)%total
!      x_i = allpoints(nJ)%center(i,1)
!      z_i = allpoints(nJ)%center(i,2)
!      e = allpoints(nJ)%layer(i)
!      
!      p = sin(Theta * PI / 180.0) / ALFA(e)
!      q = cos(Theta * PI / 180.0) / ALFA(e)
!      k = cOME * p ! this wavenumber has nothing to do with DWN summ
!      
!      ! en fortran la matriz se llena por columnas:
!      subMatSB = & 
!      RESHAPE((/ -LAMBDA(e)*p, AMU(e)*q, -p*(LAMBDA(e)+2.0*AMU(e)),&
!              -q*(LAMBDA(e)+2.0*AMU(e)), AMU(e)*p, -q*LAMBDA(e) /),&
!                          (/3,2/)) 
!      subMatSB(:,1) = subMatSB(:,1)*sin(Theta * PI / 180.0)
!      subMatSB(:,2) = subMatSB(:,2)*cos(Theta * PI / 180.0)
!      
!      ! w displacement
!      allpoints(nJ)%incident(1,i)= Uo(nJ) * & 
!                    (-1.0 * cos(Theta * PI /180.0))* &
!                    exp(-UI*k*x_i)* & 
!                    exp(UI*cOME*q*(z_i - Zstart))
!      
!      ! u displacement
!      allpoints(nJ)%incident(2,i)= Uo(nJ) * sin(Theta * PI / 180.0)* &
!                    exp(-UI*k*x_i)* & 
!                    exp(UI*cOME*q*(z_i - Zstart))
!      
!      !sigma zz
!      allpoints(nJ)%incident(3,i)= (Uo(nJ)*exp(-UI*k*x_i)* &
!                    exp(UI*cOME*q*(z_i - Zstart))*UI*cOME)* &
!                    sum(subMatSB(1,:))
!      
!      !sigma zx
!      allpoints(nJ)%incident(4,i)= (Uo(nJ)*exp(-UI*k*x_i)* &
!                    exp(UI*cOME*q*(z_i - Zstart))*UI*cOME)* &
!                    sum(subMatSB(2,:))
!      
!      !sigma xx
!      allpoints(nJ)%incident(5,i)= (Uo(nJ)*exp(-UI*k*x_i)* &
!                    exp(UI*cOME*q*(z_i - Zstart))*UI*cOME)* &
!                    sum(subMatSB(3,:))
!     end do
      
      
      subroutine matrixA_borderCond(this_A,k,cOME_i,outpf)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : verbose,UI,UR,PI!,dp
!     use waveVars, only : Theta     
!     use refSolMatrixVars, only : A ! aquí guardamos el resultado 
      use debugStuff  
!     use waveNumVars, only : K !DWN compliant / explicitly
      implicit none
      
!     real, parameter    :: x_i = 0
      complex*16,    intent(inout), dimension(4*N+2,4*N+2) :: this_A
      real               :: z_i
      real,       intent(in)     :: k
      complex*16, intent(in)     :: cOME_i  
      
      integer, intent(in) :: outpf

      complex*16, dimension(2,4) :: subMatD, subMatS
      complex*16, dimension(4,1) :: diagMat
      complex*16 :: gamma,nu,xi!,eta
      complex*16 :: egammaN,enuN,egammaP,enuP
      integer    :: iR,iC,e,bord
      
          iR= 0
          iC= 0 
          this_A = cmplx(0.0,0.0,8)
      !la matriz global se llena por columnas con submatrices de 
      !cada estrato
      DO e = 1,N+1
          
          ! algunas valores constantes para todo el estrato
          gamma = sqrt(cOME_i**2/ALFA(e)**2 - k**2)
          nu = sqrt(cOME_i**2/BETA(e)**2 - k**2)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z es mayor que cero.
          if(aimag(gamma).gt.0.0)gamma = -gamma
          if(aimag(nu).gt.0.0)nu=-nu
          
          xi = k**2 - nu**2
!         eta = 2*gamma**2 - cOME_i**2/BETA(e)**2
          
          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1 
        do bord = 0,1
          if (e+bord > N+1) then ! si 1+0;1+1;2+0;[2+1] > 2
            exit
          end if
          subMatD = RESHAPE((/ -gamma,-k*UR,-k*UR,nu,& 
                                gamma,-k*UR,-k*UR,-nu /), &
                           (/ 2,4 /))
          subMatS = RESHAPE((/ xi,-2*k*gamma,-2*k*nu,-xi,& 
                               xi,2*k*gamma,2*k*nu,-xi /),&
                           (/2,4/))  
                           
          ! la profundidad z de la frontera superior del estrato
                z_i = Z(e+bord)   ! e=1 , bord=0  ->  z = z0 = 0
                                  ! e=1 , bord=1  ->  z = Z1 = h1
          
          !downward waves
          egammaN = exp(-UI * gamma * (z_i-Z(e)))
          enuN = exp(-UI * nu * (z_i-Z(e)))
          !upward waves    
          if (e /= N+1) then !(radiation condition)
            egammaP = exp(UI * gamma * (z_i-Z(e+1)))
            enuP = exp(UI * nu * (z_i-Z(e+1)))
          else
            egammaP = 0.0d0
            enuP = 0.0d0
          end if
          
            !la matrix diagonal
            diagMat = RESHAPE((/ egammaN, enuN, egammaP, enuP /), &
                           (/ 4,1 /))
          
          ! desplazamientos estrato i (en Fortran se llena por columnas)
          subMatD = UI * subMatD !* exp(-UI * k * x_i) 
          subMatD(:,1) = subMatD(:,1) * diagMat(1,1) !Down P
          subMatD(:,2) = subMatD(:,2) * diagMat(2,1) !Down S
          subMatD(:,3) = subMatD(:,3) * diagMat(3,1) !Up   P
          subMatD(:,4) = subMatD(:,4) * diagMat(4,1) !Up   S
          
          ! esfuerzos estrato i
          subMatS = AMU(e) * subMatS !* exp(-UI*k*x_i) 
          subMatS(:,1) = subMatS(:,1) * diagMat(1,1) !Down P
          subMatS(:,2) = subMatS(:,2) * diagMat(2,1) !Down S
          subMatS(:,3) = subMatS(:,3) * diagMat(3,1) !Up   P
          subMatS(:,4) = subMatS(:,4) * diagMat(4,1) !Up   S
          
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0) then 
           if (e /= N+1) then !(radiation condition)
            if (e/=1) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+4 ) = -subMatD
            end if
             this_A( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatS
           else
             this_A( iR-1 : iR   , iC+1 : iC+2 ) = -subMatD(:,1:2)
             this_A( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatS(:,1:2)
           end if
          end if
          
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
            this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatD
            this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatS
          end if
          
        end do !bord loop del borde inferior o superior
          iR= iR+4 
          iC= iC+4
          
      END DO !{e} loop de las macro columnas para cada estrato
          
          if (verbose >= 3) then
             call showMNmatrixZ(4*N+2,4*N+2,this_A,"A    ",outpf) 
          end if

      end subroutine matrixA_borderCond
      
      
      subroutine vectorB_force(this_B,z_f,e,fisInterf,direction,cOME,k)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO
      use gloVars, only : verbose,UR,UI,PI  
      use debugStuff   
      use resultVars, only : allpoints,NPts
      implicit none
      
      complex*16,    intent(inout), dimension(4*N+2) :: this_B
      integer,    intent(in)    :: direction,e
      real,       intent(in)    :: k,z_f
      complex*16, intent(in)    :: cOME
      logical,    intent(in)    :: fisInterf
      integer :: iIf,nInterf!,iP
      real    :: SGNz,L2M!,x_i
      complex*16 :: gamma,nu,DEN
      complex*16 :: omeAlf,omeBet!,AUX
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real,       dimension(2) :: z_loc
      complex*16, dimension(2) :: egamz,enuz
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s111,s331,s131,s113,s333,s313 
                                  !  1 para Greeni1 (horizontal),
                                  !  3 para Greeni3 (vertical)
      
      
      this_B = cmplx(0.0,0.0,8)
!     B(4*N+1,:) = (-UR ) / (2.0 * PI ); return
      
      G11=0;G31=0;G33=0
      s111=0;s331=0;s131=0;s113=0;s333=0;s313=0
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
      end if   
      
      z_loc(1) = Z(e) - z_f !downward (-)
      if (e .ne. N+1) then
        z_loc(2) = Z(e+1) - z_f !upward (+)
      else
        z_loc(2) = 0.0
      end if
      
      DEN = 4.0*PI*RHO(e)*cOME**2.0
      omeAlf = cOME**2/ALFA(e)**2.0
      omeBet = cOME**2/BETA(e)**2.0
      L2M = LAMBDA(e) + 2.0*AMU(e)
      
          ! algunas valores constantes para todo el estrato          
          gamma = sqrt(omeAlf - k**2.0)
          nu = sqrt(omeBet - k**2.0)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z es mayor que cero.
          if(aimag(gamma).gt.0.0)gamma = -gamma
          if(aimag(nu).gt.0.0)nu=-nu
          
      ! en cada interfaz (1 arriba) y (2 abajo)
      do iIf = 1,2
!         print*,"[",x,",",z_loc(iIf),"]"
          if (z_loc(iIf) .ne. 0.) then
            SGNz = z_loc(iIf) / ABS(z_loc(iIf))
          else
            SGNz = 0
          end if
          egamz(iIf) = exp(-UI*gamma*ABS(z_loc(iIf)))
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))
      
      G11(iIf) = -UI/DEN * (k**2.0/gamma*egamz(iIf)+nu*enuz(iIf)) 
      G31(iIf) = -UI/DEN * SGNz*k*(egamz(iIf)-enuz(iIf)) 
      G33(iIf) = -UI/DEN * (gamma*egamz(iIf)+k**2.0/nu*enuz(iIf))
      
      s111(iIf) = -UR/DEN * ( & 
                      (k*gamma*lambda(e)+L2M*k**3.0/gamma)* egamz(iIf)&
                    + (2.0*amu(e)*k*nu) * enuz(iIf) & 
                      ) !
      
      s331(iIf) = -UR/DEN * ( &
                    (k*gamma*L2M + lambda(e)*k**3.0/gamma)* egamz(iIf)&
                  + (-2.0*amu(e)*k*nu)* enuz(iIf) &
                    ) !
                    
      s131(iIf) = -UR/DEN * amu(e)*SGNz * ( &             
                    (2.0*k**2.0)* egamz(iIf) &
                    + (nu**2.0-k**2.0)* enuz(iIf) &
                    ) !
                    
      s113(iIf) = -UR/DEN * SGNz * ( &              
                    (k**2.0*L2M + gamma**2.0*lambda(e))* egamz(iIf) &
                  + (-2.0*amu(e)*k**2.0)* enuz(iIf) &
                    ) !
                    
      s333(iIf) = -UR/DEN * SGNz * ( &              
                    (gamma**2.0*L2M + k**2.0*lambda(e))* egamz(iIf) &
                  + (2.0*amu(e)*k**2.0)* enuz(iIf) &
                    ) !              
                    
      s313(iIf) = -UR/DEN * amu(e) * ( &              
                    (2.0*k*gamma)* egamz(iIf) &
                  - (k/nu*(nu**2.0-k**2.0))* enuz(iIf) &
                    ) !
      
      end do !iIf interface
          
      if (fisInterf) then
      if(direction .eq. 1) then
        s331(1) = (-UR ) / (2.0 * PI )
      elseif (direction .eq. 3) then
        s333(1) = (-UR ) / (2.0 * PI )
      end if
        G33(1) = 0
        G31(1) = 0
      end if
      
      ! El vector de términos independientes genera el campo difractado
      
      if (direction .eq. 1) then ! fuerza HORIZONTAL
       if (e .ne. 1) then ! =     (1) interfaz de arriba; (2) abajo
        this_B(1+4*(e-1)-2) =  G31(1)!  w
        this_B(1+4*(e-1)-1) =  G11(1)!  u
       end if 
        this_B(1+4*(e-1)  ) =  S331(1)! szz
        this_B(1+4*(e-1)+1) =  S131(1)! szx
    
      if (.not. fisInterf) then  
       if (e .ne. N+1) then! =      (1) interfaz de arriba; (2) abajo
        this_B(1+4*(e-1)+2) = - G31(2)!  w
        this_B(1+4*(e-1)+3) = - G11(2)!  u
        this_B(1+4*(e-1)+4) = - S331(2)! szz
        this_B(1+4*(e-1)+5) = - S131(2)! szx
       end if
      end if
      elseif (direction .eq. 2) then ! fuerza VERTICAL
       if (e .ne. 1) then! =     (1) interfaz de arriba; (2) abajo
        this_B(1+4*(e-1)-2) =  G33(1)!  w 
        this_B(1+4*(e-1)-1) =  G31(1)!  u 
       end if 
        this_B(1+4*(e-1)  ) =  S333(1)! szz
        this_B(1+4*(e-1)+1) =  S313(1)! szx 
    
      if (.not. fisInterf) then  
       if (e .ne. N+1) then!=      (1) interfaz de arriba; (2) abajo
        this_B(1+4*(e-1)+2) = - G33(2)!  w 
        this_B(1+4*(e-1)+3) = - G31(2)!  u
        this_B(1+4*(e-1)+4) = - S333(2)! szz 
        this_B(1+4*(e-1)+5) = - S313(2)! szx 
       end if
      end if
      
      end if ! direction
      end subroutine vectorB_force
      function calcMecaAt_xi_zi(thisIP_B,z_i,e,cOME_i,k,outpf)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : verbose,UI,UR,PI
      use waveVars, only : Theta
      use resultVars, only : MecaElem
      implicit none
      type (MecaElem)              :: calcMecaAt_xi_zi
      real, intent(in)             ::  z_i, k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e,outpf
      complex*16, dimension(4*N+2),intent(in) :: thisIP_B
      complex*16 :: gamma,nu,eta,xi
      complex*16 :: egammaN,enuN,egammaP,enuP
      complex*16, dimension(2,4) :: subMatD
      complex*16, dimension(3,4) :: subMatS
      complex*16, dimension(4,1) :: diagMat
      complex*16, dimension(4,1) :: partOfXX
      complex*16, dimension(2,1) :: resD
      complex*16, dimension(3,1) :: resS
      
      if (verbose >= 3) then
       write(outpf,'(a,F7.3,a,F12.7,a,F10.2,a,F10.2,a,I0)') & 
                    "finding solution values at w:", & 
                    real(cOME_i),"k=",k," {",z_i,"} e=",e
      end if
      
!     calcMecaAt_xi_zi%center%X = x_i
!     calcMecaAt_xi_zi%center%Z = z_i
       
      ! algunas valores constantes para todo el estrato
      gamma = sqrt(cOME_i**2. /ALFA(e)**2. - k**2.)
      nu = sqrt(cOME_i**2. /BETA(e)**2. - k**2.)
      
      if(aimag(gamma).gt.0.0)gamma = -gamma
      if(aimag(nu).gt.0.0)nu=-nu
      
      xi = k**2. - nu**2.
      eta = 2.0*gamma**2.0 - (cOME_i/BETA(e))**2.0
      subMatD = RESHAPE((/ -gamma,-k*UR,-k*UR,nu,gamma,-k*UR,-k*UR,-nu /), &
                           (/ 2,4 /))
      subMatS = RESHAPE((/ xi,      -2.0*k*gamma,     eta,     &
                          -2.0*k*nu,     -xi,        2.0*k*nu,   &
                           xi,       2.0*k*gamma,     eta,     &
                           2.0*k*nu,     -xi,       -2.0*k*nu /),&
                           (/3,4/))      
      
         !downward waves
          egammaN = exp(-UI * gamma * (z_i-Z(e)))
          enuN = exp(-UI * nu * (z_i-Z(e)))
          !upward waves 
          if (e /= N+1) then !(radiation condition)
            egammaP = exp(UI * gamma * (z_i-Z(e+1)))
            enuP = exp(UI * nu * (z_i-Z(e+1)))
          else
            egammaP = 0.0d0
            enuP = 0.0d0
          end if
          !la matrix diagonal
            diagMat = RESHAPE((/ egammaN, enuN, egammaP, enuP /), &
                           (/ 4,1 /))
          ! desplazamientos estrato i (en Fortran se llena por columnas)
!         subMatD = UI * exp(-UI * k * x_i) * subMatD
          subMatD = UI * subMatD
          subMatD(:,1) = subMatD(:,1) * diagMat(1,1) !Down P
          subMatD(:,2) = subMatD(:,2) * diagMat(2,1) !Down S
          subMatD(:,3) = subMatD(:,3) * diagMat(3,1) !Up   P
          subMatD(:,4) = subMatD(:,4) * diagMat(4,1) !Up   S
          
          ! esfuerzos estrato i
!         subMatS = AMU(e) * exp(-UI * k * x_i) * subMatS
          subMatS = AMU(e) * subMatS
          subMatS(:,1) = subMatS(:,1) * diagMat(1,1) !Down P
          subMatS(:,2) = subMatS(:,2) * diagMat(2,1) !Down S
          subMatS(:,3) = subMatS(:,3) * diagMat(3,1) !Up   P
          subMatS(:,4) = subMatS(:,4) * diagMat(4,1) !Up   S
      
 !mutiplicamos por los coeficientes de la solución del sistema
        if (e /= N+1) then
          partOfXX(1:4,1) = thisIP_B(4*(e-1)+1 : 4*(e-1)+4)
        else !( condición de radiación)
          partOfXX(1:2,1) = thisIP_B(4*(e-1)+1 : 4*(e-1)+2)
          partOfXX(3:4,1) = 0.0d0
        end if!
        if (verbose>=3) then
          print*,"xx1 ",partOfXX(1,1)
          print*,"xx2 ",partOfXX(2,1)
          print*,"xx3 ",partOfXX(3,1)
          print*,"xx4 ",partOfXX(4,1)
        end if
        resD = matmul(subMatD,partOfXX)
        resS = matmul(subMatS,partOfXX)
        calcMecaAt_xi_zi%Rw(1) = resD(1,1)
        calcMecaAt_xi_zi%Rw(2) = resD(2,1)
        calcMecaAt_xi_zi%Rw(3) = resS(1,1)
        calcMecaAt_xi_zi%Rw(4) = resS(2,1)
        calcMecaAt_xi_zi%Rw(5) = resS(3,1)
          
      end function calcMecaAt_xi_zi
      
      !     subroutine Add_diffractedField_by_layeredMedia_source & 
!     (NumObservers,P,theseIPs_B,direction,cOME_i,iik,outpf)
!     use gloVars, only : imecMax
!     use resultVars, only : MecaElem,Punto!,allpoints,NPts
!     use soilVars, only : N
!     
!     implicit none
!     integer, intent(in)          :: NumObservers,iik,outpf
!     integer,    intent(in)       :: direction
!     complex*16, intent(in)       :: cOME_i
!     type(Punto), dimension(NumObservers), intent(inout) :: P
!     complex*16, dimension(4*N+2),intent(in) :: theseIPs_B
!     type(MecaElem)               :: calcMecaAt_xi_zi, this_Meca
!     real                         :: x_i, z_i
!     integer                      :: iP,e
!     do iP = 1,NumObservers
!       x_i = P(iP)%center(1)
!       z_i = P(iP)%center(2)
!         e = P(iP)%layer
!       this_Meca = calcMecaAt_xi_zi(theseIPs_B,x_i,z_i,e,cOME_i,outpf)
!       
!       if (direction .eq. 1) then !x
!         P(iP)%FKh(iik,1:imecMax) = & 
!         P(iP)%FKh(iik,1:imecMax) + this_Meca%Rw(1:imecMax)
!       elseif (direction .eq. 2) then !z
!         P(iP)%FKv(iik,1:imecMax) = & 
!         P(iP)%FKv(iik,1:imecMax) + this_Meca%Rw(1:imecMax)
!       end if
 !     end do
!     end subroutine Add_diffractedField_by_layeredMedia_source
!     
!     
!     subroutine Add_diffractedField_by_layeredMedia_Green & 
!     (NumObservers,P,theseIPs_B,iPxi,iGq,direction,cOME_i,iik,outpf)
!     use gloVars, only : imecMax
!     use resultVars, only : MecaElem,Punto!,allpoints,NPts
!     use soilVars, only : N
!     
!     implicit none
!     integer, intent(in)          :: NumObservers,iPxi,iGq,direction,iik,outpf
!     complex*16, intent(in)       :: cOME_i
!     type(Punto), dimension(NumObservers), intent(inout) :: P
!     complex*16, dimension(4*N+2),intent(in) :: theseIPs_B
!     type(MecaElem)               :: calcMecaAt_xi_zi, this_Meca
!     real                         :: x_i, z_i
!     integer                      :: iP,e
!     do iP = 1,NumObservers
!       x_i = P(iP)%center(1)
!       z_i = P(iP)%center(2)
!         e = P(iP)%layer
!       this_Meca = calcMecaAt_xi_zi(theseIPs_B,x_i,z_i,e,cOME_i,outpf)
!       
!       !guardar en: GT_gq_k(xi,iGq,imec,m,k)
!       P(iP)%GT_gq_k(iPxi,iGq,1:imecMax,direction,iik) = &
!       P(iP)%GT_gq_k(iPxi,iGq,1:imecMax,direction,iik) + this_Meca%Rw(1:imecMax)
!       !                                 1 == horizontal
!       !                                 2 == vertical
!     end do 
!     
!     end subroutine Add_diffractedField_by_layeredMedia_Green      
      
      
      function calcFF(z_f,e_f,dir,zX,eX,k,cOME)
      use soilVars 
      use gloVars  
      use waveNumVars, only : NFREC,NMAX,DK
      use resultvars, only: Punto,MecaElem
      
      implicit none
      integer,    intent(in)     :: e_f,eX,dir
      real,       intent(in)     :: z_f,zX
      real,       intent(in)     :: k
      complex*16, intent(in)     :: cOME
      type (MecaElem)            :: calcFF
      
      real                       :: z_i,SGNz
      complex*16 :: G11,G33,G31,S333,S313,S113,s331,s131,s111
      complex*16 :: egamz,enuz
      real                       :: L2M
      complex*16                 :: gamma,nu,DEN
      complex*16                 :: omeAlf,omeBet
      
      egamz=cmplx(1.0,0.0,8);enuz=cmplx(1.0,0.0,8);SGNz=0
      calcFF%Rw = cmplx(0.0,0.0,8)
      
      if (e_f .ne. eX) return
      
      DEN = 4.0*PI*RHO(e_f)*cOME**2.0
      omeAlf = cOME**2/ALFA(e_f)**2.0
      omeBet = cOME**2/BETA(e_f)**2.0
      L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
      
      z_i = zX - z_f
      if (z_i .ne. 0) then
        SGNz = z_i / ABS(z_i)
      end if
      
        gamma = sqrt(omeAlf - k**2.0)
        nu = sqrt(omeBet - k**2.0)
        if(aimag(gamma).gt.0.0)gamma = -gamma
        if(aimag(nu).gt.0.0)nu=-nu
         
          egamz = exp(-UI*gamma*ABS(z_i))
          enuz = exp(-UI*nu*ABS(z_i)) 
          
          
          G31 = -UI/DEN * SGNz*k*(egamz - enuz)
          
      if(dir .eq. 2) then
         G33 = -UI/DEN * (gamma*egamz + k**2.0/nu*enuz) 
         S333 = -UR/DEN * SGNz*( &
                    (gamma**2.0*L2M + k**2.0*lambda(e_f))* egamz &
                  + (2.0*amu(e_f)*k**2.0)* enuz &
                    ) 
         S313 = -UR/DEN * amu(e_f) * ( &
                    (2.0*k*gamma)* egamz &
                  - (k/nu*(nu**2.0-k**2.0))* enuz &
                    ) 
         S113 = -UR/DEN * SGNz * ( &
                    (k**2.0*L2M + gamma**2.0*lambda(e_f))* egamz &
                  + (-2.0*amu(e_f)*k**2.0)* enuz &
                    ) 
      
        calcFF%RW(1) = G33
        calcFF%RW(2) = G31
        calcFF%RW(3) = S333
        calcFF%RW(4) = S313
        calcFF%RW(5) = S113
      else
         G11 = -UI/DEN * (k**2.0/gamma*egamz + nu*enuz) 
         
         s331 = -UR/DEN * ( &
                    (k*gamma*L2M + lambda(e_f)*k**3.0/gamma)* egamz &
                  + (-2.0*amu(e_f)*k*nu)* enuz &
                    )
         s131 = -UR/DEN * amu(e_f)*SGNz * ( &             
                    (2.0*k**2.0)* egamz &
                    + (nu**2.0-k**2.0)* enuz &
                    )
         s111 = -UR/DEN * ( & 
                    (k*gamma*lambda(e_f)+L2M*k**3.0/gamma)* egamz &
                    + (2.0*amu(e_f)*k*nu)* enuz & 
                    )
     
        calcFF%RW(1) = G31
        calcFF%RW(2) = G11
        calcFF%RW(3) = S331
        calcFF%RW(4) = S131
        calcFF%RW(5) = S111
      end if
!       calcFF%RW = calcFF%RW * exp(-UI * k * xX)
      end function calcFF            
      subroutine inverseA(A,n)
      integer, intent(in) :: n
      complex*16, dimension(n,n), intent(inout) :: A
      
      integer, dimension(:),allocatable :: ipiv
      integer :: lwork
      complex*16, dimension(:),allocatable :: work
      integer :: info
      
      allocate(ipiv(n+1))
      lwork = n*n
      allocate(work(lwork))
      call zgetrf(n,n,A,n,ipiv,info)
      if(info .ne. 0) stop "Problem at LU factorization of matrix "
      call zgetri(n,a,n,ipiv,work,lwork,info)
      if(info .ne. 0) stop "Problem at inverse of matrix "
      
      deallocate(work)
      deallocate(ipiv)
      end subroutine inverseA
      subroutine makeTaperFuncs(Npol,cutoff_fracFmax,cutoff_fracKmax)
      use resultvars, only: Hf,Hk 
      use waveNumVars, only : NFREC,NMAX,DFREC,DK
      
      implicit none
      integer, intent(in) :: Npol
      real, intent(in) :: cutoff_fracFmax, cutoff_fracKmax
      integer :: i
      
      allocate(Hf(NFREC+1))
      allocate(Hk(NMAX*2))
      
         Hf = (/ (i,i=0,nfrec) /) * DFREC
         Hf = 1.0/sqrt(1.0+(Hf/(cutoff_fracFmax*NFREC*DFREC))**(2*Npol)) 
         
         Hk(1:nmax) = (/ (i,i=0,nmax-1) /) * DK
         Hk(nmax+1:2*nmax) = (/ (nmax-i,i=0,nmax-1) /) * DK 
         HK = 1.0/sqrt(1.0+(Hk/(cutoff_fracKmax*NMAX*DK))**(2*Npol))
      end subroutine makeTaperFuncs  
      
      subroutine expPIK(expK)
      use waveNumVars, only : NMAX,DK
      use glovars, only : PI
      use debugStuff
      implicit none
      complex*16, dimension(2*nmax) :: expK
      real*8, dimension(2*nmax) :: aux
      integer :: i
      
      aux(1:nmax) = (/ (i,i=0,nmax-1) /) * DK
      aux(1) = dk * 0.001
      aux(nmax+1:2*nmax) = (/ (nmax-i,i=0,nmax-1) /) * DK * (-1)
      expK = exp(cmplx(0.0 ,-pi * aux,8))
      
!     call showMNmatrixZ (2*nmax,1,expK,"nmax ",6)
      end subroutine expPIK
      function freefTraction(iP_x,l)
      use resultVars, only:BouPoints!,nBpts
      implicit none
      complex*16 :: freefTraction
      integer, intent(in) :: iP_x,l
      
      !       en realidad ya es FX ---, ,--- nomas hay una, la deorita
      !                               | |
      freefTraction = BouPoints(iP_x)%W(1,5-l) * BouPoints(iP_x)%normal(1) + &
                      BouPoints(iP_x)%W(1,4-l) * BouPoints(iP_x)%normal(2)
!     print*,BouPoints(iP_x)%normal(1)," ",BouPoints(iP_x)%normal(2)
      ! las tracciones:
      !
      !   ,--- componente de la tracción : m
      !   |     ,--- cara
      !   |     |
      !  Tx = Sxx nx + Szx nz  |___ (0) fza real (fuente)  .
      !  Tz = Szx nx + Szz nz  |                           .  
      !
      !  T_lm = S_lkm * n_k
      !       = S_l1m * n_1 + S_l3m * n3
      !         __|__         __|__
      !      s11     s31   s13    s33
      !      s11     s31   s31    s33  (son equivalentes)
      !       5       4     4      3   (indice en w(_,i) )
      !       0       1     0      1   (indice de submatriz: l )
      end function freefTraction
      function integralEq16(iP_x,l,iPxi,m)
         !                              receptor---,       ,--- fuente
         !                                        _|___   _|___
         !  ibemMat(iP_x+l,iPxi+m) = integralEq14(iP_x,l,iPxi,m)
      use resultVars, only: P => allpoints, B => boupoints
                     
      use Gquadrature, only: Gquad_n
      implicit none
      complex*16 :: integralEq16,desp
      integer, intent(in) :: iP_x,l,iPxi,m
      logical :: lejos
      integer :: iGq
      
      integralEq16 = 0
      lejos = .false.
      
      if (lejos) then
      ! no usamos integración gaussiana
             
      else ! cerca
      ! usamos integracion gaussiana
      
       ! En cada punto gaussiano del segmento en Xi en la dirección l
       do iGq = 1,Gquad_n
!        print*,"iGq=",iGq
         ! dada la fuerza en XI en la dirección (m+1)
         ! desplazamiento en la dirección (l+1)
         !     receptor        emisor    
         desp = P(iP_x)%GT_gq(iPxi,iGq,2-l,m+1)  
!        print*,"desp=",desp
         ! multiplicar por peso gaussiano de integración en segmento xi 
         desp = desp * B(iPxi)%Gq_xXx_C(iGq)
!        print*,"peso=",B(iPxi)%Gq_xXx_C(iGq)
         ! acumular suma
         integralEq16 = integralEq16 + desp
       end do !xXx
      end if !lejos
!     print*,integralEq16 
      end function integralEq16
      
      ! los deplazamientos
      ! iMec     U        W
      !  i       2        1
      !  l       0        1
      !necesito  x        z
      function integralEq16mov(iP,iP_x,l,iPxi,m)
         !                              receptor---,       ,--- fuente
         !                                        _|___   _|___
         !  ibemMat(iP_x+l,iPxi+m) = integralEq14(iP_x,l,iPxi,m)
      use resultVars, only:P => allpoints, B => boupoints
      use Gquadrature, only: Gquad_n
      implicit none
      complex*16 :: integralEq16mov,desp
      integer, intent(in) :: iP,iP_x,l,iPxi,m
      logical :: lejos
      integer :: iGq
      
      integralEq16mov = 0
      lejos = .false.
      
      if (lejos) then
      ! no usamos integración gaussiana
             
      else ! cerca
      ! usamos integracion gaussiana
      
       ! En cada punto gaussiano del segmento en Xi en la dirección l
       do iGq = 1,Gquad_n
         ! dada la fuerza en XI en la dirección (m+1)
         ! desplazamiento en la dirección (l+1)
         !     receptor        emisor    
         desp = P(iP_x)%GT_gq_mov(iPxi,iGq,2-l,m+1,iP) 
         ! multiplicar por peso gaussiano de integración en segmento xi 
         desp = desp * B(iPxi)%Gq_xXx_C(iGq)
         ! acumular suma
         integralEq16mov = integralEq16mov + desp
       end do !xXx
      end if !lejos
      
      end function integralEq16mov
      
      ! los deplazamientos
      ! iMec     U        W
      !  i       2        1
      !  l       0        1
      !necesito  x        z
      subroutine W_to_t(W,iP,coords,outpf)
      use waveNumVars, only : NFREC,DFREC
      use glovars
      use waveVars, only : dt,t0,Uo
      use ploteo10pesos
      use wavelets
      implicit none
      integer ,intent(in) :: iP,outpf
      complex*16, dimension(Nfrec+1,2), intent(in) :: W
      complex*16, dimension(2*NFREC) :: S
      real, dimension(2), intent(in) :: coords  
!     character(LEN=10),intent(in) :: tt
      character(LEN=3), dimension(5)   :: nombre
      character(LEN=100) :: titleN
      character(LEN=9)   :: logflag
      integer :: i,iMec,t
      real :: x_i,z_i,factor
      factor = sqrt(2.*NFREC)
      
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 's33'
      nombre(4)= 's31'
      nombre(5)= 's11'
      
      x_i=coords(1)
      z_i=coords(2)
      
      do iMec = 1,2
      !  (0) crepa en f
      !           1 2 3 4 ... Nfrec Nfrec+1
      ! tenemos:  0 1 2 3 ... N/2-1  N/2
      ! hacemos:  0 1 2 3 ... N/2-1 -N/2 ... -3 -2 -1
      S(      1:nfrec  ) =       W(1:nfrec:+1,  iMec) * & 
      exp(cmplx(0.0,- pi * (/((t-1)*dfrec,t=1,nfrec)/) * t0,8)) * Uo(1:nfrec)
      S(nfrec+1:nfrec*2) = conjg(W(nfrec+1:2:-1,iMec)) * & 
      exp(cmplx(0.0,- pi * (/((-t)*dfrec,t=nfrec,1,-1)/) * t0,8)) * Uo(nfrec+1:nfrec*2)
      
      
      if (verbose .ge. 2) then
      !guardar en texto
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
!        print*,titleN
         OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         write(3211,'(I0)') 1
         write(3211,'(a)') "amplitude"
         write(3211,'(F15.8)') DFREC
         do i = 1,size(S)
          write(3211,'(ES14.5E2,2x,ES14.5E2)') real(S(i)),aimag(S(i))
         end do
         close (3211) 
      end if
      
      ! imprimir espectro:
        if (Verbose .ge. 2) then
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
!        print*,titleN
         call plotXYcomp(S,DFREC,2*nfrec,titleN, & 
         'frec[hz] ','amplitude',1200,800)
         
        end if
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'fL_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         logflag = 'logx     '
         call plotSpectrum(S,DFREC,2*nfrec,nfrec, & 
                         titleN,'frec[hz]','amplitude',logflag,1200,800)

      !  (1) pasar al tiempo
         S = S*factor
         call fork(2*nfrec,S,+1,verbose,outpf)
         S = S/factor
      !  (2) remover efecto de la frecuencia imaginaria
         S = S * exp(-DFREC/periodicdamper * dt*((/(i,i=1,2*nfrec)/)-1))
      
      !  (3) plot the damn thing
      if (verbose .ge. 1) then
      !guardar en texto
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'S_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
!        print*,titleN
         OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         write(3211,'(I0)') 1
         write(3211,'(a)') "amplitude"
         write(3211,'(F15.8)') DT
         do i = 1,size(S)
          write(3211,'(ES14.5E2,2x,ES14.5E2)') real(S(i)),aimag(S(i))
         end do
         close (3211) 
      end if
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'S_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S,dt,2*nfrec,titleN, & 
         'time[sec]','amplitude',1200,800)
      
      end do !imec
      end subroutine W_to_t
      subroutine Hollywood(outpf)
      use DISLIN
      use resultVars
      use wavelets ! FORK(LX,CX,SIGNI,verbose,outpf)
      use waveVars, only : Dt,Uo,t0
      use waveNumVars ! FK,NFREC,DFREC
      use refSolMatrixVars, only : COME
      use gloVars
      use MeshVars
      use soilVars, only : Z,N
      use GeometryVars, only : Xcoord,nxi
      implicit none
      integer ,intent(in) :: outpf
      real       ::  X(nxXxmpts)
      real       ::           Y(nz)
      complex*16 :: Sm(nxXxmpts,nz,imecMax,2*NFREC)      
      integer  :: iP,iMec
      integer  :: ix,iz,iT,i,Iframe,Fframe,mecMAx,t
      real     :: factor,ColorRangeMaximumScale,tlabel
      integer*4 :: lentitle
      character(LEN=3), dimension(5)      :: nombre
      character(LEN=60) :: textoTimeStamp
      character(LEN=15) :: auxTxt
      CHARACTER(len=400) :: path
      character(LEN=1) :: imdone
      character(LEN=60) :: CTIT

      real :: maxY, minY, maxX, minX, p
      real, dimension(41)   :: ZLVRAY
      real                  :: maV,miV,Vstep,xstep,zstep
      
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 's33'
      nombre(4)= 's31'
      nombre(5)= 's11'
      mecMax = 2
      if(workboundary) mecMax = 2
      
      if (verbose >= 1) Write(outpf,'(a)') "Will make a movie..."
         CALL chdir("video")
         CALL getcwd(path)
         write(outpf,'(a,a)') "at ",TRIM(path)
         call system("rm *.png")
         call system("rm *.txt")
      
      factor = sqrt(real(2.*NFREC))
      
      allocate(Sismogram(size(X),size(Y),5,2*NFREC))
      Sismogram = 0
      ! (0) salir de K --- listo en FKtoFX
      ! (1) reordenar creapa exitente en X para encontrarlos más facil
!        FKm = cshift(FKm,SHIFT=nmax/2+1,DIM=1)
      ! (1.1) coordenadas X 
      i = 1
!     print*,size(x)
      do ix=-MeshDXmultiplo * (NxXxMpts-1)/2, MeshDXmultiplo * (NxXxMpts-1)/2, MeshDXmultiplo
        X(i) = ix * DeltaX  !;print*,'x=',X(i)
        i = i + 1
      end do
      
      ! se hace para los puntos del video.
      iz = 1
      do iP= mPtini,mPtfin 
      do ix = 1,NxXxMpts
        Y(iz) = allpoints(IP)%center(2) !;print*,iP,'z=',Y(iz)
        
        ! los espectros de cada coordenada ya están ordenaditos
        do iMec = 1,mecMax
        ! (1) crepa
        Sm(ix,iz,iMec,1:nfrec) = allpoints(iP)%WmovieSiblings(ix,1:nfrec:+1,iMec)* & 
      exp(cmplx(0.0,- pi * (/((t-1)*dfrec,t=1,nfrec)/) * t0,8)) * Uo(1:nfrec)
      
        Sm(ix,iz,iMec,nfrec+1:nfrec*2) = conjg(allpoints(iP)%WmovieSiblings(ix,nfrec+1:2:-1,iMec))* & 
      exp(cmplx(0.0,- pi * (/((-t)*dfrec,t=nfrec,1,-1)/) * t0,8)) * Uo(nfrec+1:nfrec*2)
        
        ! (2) pasar al tiempo: fork
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) * factor
          call fork(2*nfrec,Sm(ix,iz,iMec,:),+1,verbose,outpf)
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) / factor
        
        ! (2.3) remover efecto de la frecuencia imaginaria
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) * & 
                            exp(-DFREC/periodicdamper * Dt*((/(i,i=1,2*nfrec)/)-1))
        
        end do !iMec
      end do !ix
      iz = iz + 1
      end do !iP
      
      !color table boundaries
      ColorRangeMaximumScale = 0.1
      
  123 do i=1,mecMax
       maV = maxVal(real(Sm(:,:,i,:)))
       miV = minVal(real(Sm(:,:,i,:)))
       maV = max(maV,abs(miV))
       miV = - maV
       
       print *, char(7)
       write(6,'(a,a,a,F10.6,a,/,a,/,a,F10.6,/,a)', ADVANCE = "NO") & 
      'Look at the seismograms for ', nombre(i), & 
      '. Is the response too spiky (the max = ',maV,'? ', &
      'We can enhance detail by reducing the value for maximum color. ', &
      'We propose the maximum to be = ', & 
      ColorRangeMaximumScale * maV, &
      'Do yo want to change it [Y]/[N] ?'
      
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        write(6,'(a)', ADVANCE = "NO")'New maximum for plot = '
        read(5,*) p
        ColorRangeMaximumScale = p / maV
      end if
      
       colorBounds(i,1) = maV * ColorRangeMaximumScale
       colorBounds(i,2) = miV * ColorRangeMaximumScale
       
       if (verbose >= 1) then
        write(outpf,'(a,a,a,a,E12.4,a,E12.4,a)') "colorbounds:", & 
        "(",nombre(i),") [",colorBounds(i,2)," ; ",colorBounds(i,1),"]"
       end if
      end do
!     stop 0
            ! plotting boundaries according to geometry of topography
      minX=0.;maxX=0.;minY=0.;maxY=0.
      if (workBoundary) then
         minX = MINVAL(Xcoord(:,1))
         maxX = MAXVAL(Xcoord(:,1))
         minY = MINVAL(Xcoord(:,2))
         maxY = MAXVAL(Xcoord(:,2))
      end if
      minX = MIN(MIN(minX,-1.0),MINVAL(X))!; print*,minX
      maxX = MAX(MAX(maxX, 1.0),MAXVAL(X))!; print*,maxX
      minY = MIN(MIN(minY, 0.0),MINVAL(Y))!; print*,minY
      maxY = MAX(MAX(maxY, 2.0),MAXVAL(Y))!; print*,maxY
      
      Iframe = 1
      Fframe = size(Sm(1,1,1,:))
      print *, char(7)
      write(6,'(a,I0,a)')'Look at the seismograms, I will plot ',Fframe,&
      'frames. Proceed [Y] or change frame range [else]'
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        Write(outpf,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
      else
        write(6,'(a)', ADVANCE = "NO")'Start frame = '
        read(5,*)Iframe
        write(6,'(a)', ADVANCE = "NO")'Final frame = '
        read(5,*)Fframe
        Write(outpf,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
      end if
      
      
      do iT=Iframe,Fframe !cada fotograma
      do iMec=1,mecMax !cada elemento mecanico 
       write(auxTxt,'(F8.3,30x)') Dt*(iT-1)
       if (verbose >= 2) then
        write(textoTimeStamp,'(a,a,a,a)') & 
         nombre(iMec), '-' ,trim(adjustl(auxTxt)),'.txt'
        open(44, file=textoTimeStamp)
        do iz=1,nz
        do ix=1,size(Sm,1)
          write(44,'(e15.6)',advance='no') real(Sm(ix,iz,iMec,iT))
        end do
        write(44,*)
        end do
        close(44)
       end if
      
!     write(textoTimeStamp,'(a,a,a)') & 
!        nombre(iMec), trim(adjustl(auxTxt)),'.pdf'
!     write(outpf,'(a)') textoTimeStamp
      write(textoTimeStamp,'(a,I0,a)') nombre(iMec), iT,'.png'
      
      maV = colorBounds(iMec,1)
      miV = colorBounds(iMec,2)
      Vstep = (maV-miV)/40.0
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      
      ! shaded contour plot
      CALL METAFL('PNG') !'PDF'
!     print*,'will print ', trim(textoTimeStamp)
      CALL SETFIL(trim(textoTimeStamp))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(3600,4),int(2400,4))
      call imgfmt('RGB')
      call winsiz(int(1200),int(800,4))
      CALL SCRMOD('REVERS') !fondo blanco
      CALL DISINI()
!     CALL COMPLX ! sets a complex font
      CALL BMPFNT ('SIMPLEX')
!     CALL HWFONT()
           !the position of an axis system. Lower left corner
      CALL axspos (int(300,4) ,int(2200,4)) 
      call axslen (int(2700,4), int(2000,4)) !size of the axis system.
      
      xstep = abs(X(1))/3.0
      
!     xstep = max(X(2)-X(1),real(int((maxX-minX) / 5.0 )))
      zstep = max(Y(2)-Y(1),real(int( (maxY-minY) / 10.)))
      
      !number of decimal places for labels
      call labdig(int(2,4),'Z') 
      if (xstep .gt. 10) then; call labdig(int(-1,4),'X')
        else; call labdig(int(1,4),'X'); end if!
      if (zstep .gt. 10) then; call labdig(int(-1,4),'Y')
        else; call labdig(int(1,4),'Y'); end if!
      
      call setgrf("TICKS", "NAME", "NAME", "TICKS")      
      CALL LABELS ('EXP ','Z')
      CALL SETVLT ('SPEC')
      
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4),& 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
      
      CALL CONSHD(real(X,4),int(size(X),4),real(Y,4),int(size(Y),4), & 
                real(real(Sm(:,:,iMec,iT)),4), real(ZLVRAY,4),int(41,4))
      
      call color ('FORE')
      call shdpat(int(0,4))
      do i=1,N
         call rlrec(real(minX,4),real(Z(i),4), & 
                    real(maxX-minX,4),real(Z(i+1)-Z(i),4))
      end do
      
      if (workboundary) then
      ! the topography
      call color ('BACK')
      call shdpat(int(16,4))
      call rlarea(real(Xcoord(:,1),4),real(Xcoord(:,2),4),int(nXI,4))
      call color ('FORE')
      
      CALL HSYMBL(int(7,4))
      CALL MRKCLR(int(-1,4))
      call marker(int(15,4))
      
      call curve(real(Xcoord(:,1),4),real(Xcoord(:,2),4),int(nXI,4))
!     call xaxgit()
      
      call color ('BLACK')
      call rline(real(Xcoord(1,1),4),real(Xcoord(1,2),4),real(Xcoord(nXI,1),4),real(Xcoord(nXI,2),4))
      call color ('FORE')
      end if
      
!     CTIT='Vectors'
      tlabel = (it)*dt
!     print*,tlabel
      write(CTIT,'(a,F9.5,a)') 't=',tlabel,' seg'
      lentitle = NLMESS(CTIT)
      CALL MESSAG(CTIT,int((3600-lentitle-100),4),int(300,4))
      
      CALL HEIGHT(int(150,4))
      call errmod ("all", "off")
      CALL DISFIN 
      end do !iMec
!     stop "killed video"
      end do !iT
      
      !now make the video with the frames
      do iMec=1,mecMax
      write(path,'(a,a,a)')'ffmpeg -i ',nombre(iMec), & 
                  '%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4'
      call system(trim(path))
      write(path,'(a,a,a)') 'cp video.mp4 ',nombre(iMec),'.mp4'
      call system(trim(path))
      call system('rm video.mp4')
      end do
      
      print *, char(7)
      write(6,'(a)') 'video is done. Do you want to change it?'
      write(6,'(a)', ADVANCE = "NO") 'replot [Y] , no I am done [else]: '  
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        go to 123
      end if
      print *, char(7)
      write(6,'(a)') 'Delete *.png files? [Y]'
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        write(path,'(a)') 'rm *.png'
        call system(trim(path))
      end if
      if (verbose .ge. 2) then
        write(6,'(a)') 'Delete *.txt files? [Y]'
        read(5,*)imdone
        if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
          write(path,'(a)') 'rm *.txt'
          call system(trim(path))
        end if
      end if
      end subroutine Hollywood

