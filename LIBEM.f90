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
      integer    :: multSubdiv ! Lo ideal es 6 
      real       :: periodicdamper
      logical    :: plotStress
      logical    :: plotModule
      logical, parameter    :: getInquirePointsSol = .true.!dont change
      integer    :: imecMax ! 2 o 4
      integer    :: plusOne ! 0 o 1
      real, dimension(2), parameter :: fixedNormal = (/0.0,1.0/)
      complex*16, parameter :: UI = cmplx(0.0,1.0,8), &
                               UR = cmplx(1.0,0.0,8), &
                               Z0 = cmplx(0.0,0.0,8)
      real, parameter :: PI = real(4.0*ATAN(1.0),8)
      integer, parameter :: Hplot = 700 , Wplot = 1200
      end module gloVars
      
      module Gquadrature
      real :: WLmulti ! una o dos longitudes de onda
      integer, parameter :: Gquad_n = 3
      ! ##############################
      ! # Para cambiar Num de puntos #
      ! # actualizar encabezado de   #
      ! # rutinas de topografía      #
      ! ##############################
      
      ! courtesy of http://pomax.github.io/bezierinfo/legendre-gauss.html
      real, parameter, dimension(3) :: Gqu_t_3 = & 
      (/ -0.7745966692414834, &
          0.0000000000000000, &
         +0.7745966692414834 /)
      
      real, parameter, dimension(3) :: Gqu_A_3 = & 
      (/  0.5555555555555556, &
          0.8888888888888888, &
          0.5555555555555556 /)
      
!     real, parameter, dimension(7) :: Gqu_t_7 = & 
!     (/ -0.9491079123427585, &
!        -0.7415311855993945, &
!        -0.4058451513773972, &
!         0.0000000000000000, &
!        +0.4058451513773972, &
!        +0.7415311855993945, &
!        +0.9491079123427585 /)
!     
!     real, parameter, dimension(7) :: Gqu_A_7 = & 
!     (/  0.1294849661688697, &
!         0.2797053914892766, &
!         0.3818300505051189, &
!         0.4179591836734694, &
!         0.3818300505051189, &
!         0.2797053914892766, &
!         0.1294849661688697 /)
      
      real, parameter, dimension(8) :: Gqu_t_8 = & 
      (/ -0.9602898564975363, &
         -0.7966664774136267, &
         -0.5255324099163290, &
         -0.1834346424956498, &
          0.1834346424956498, &
          0.5255324099163290, &
          0.7966664774136267, &
          0.9602898564975363 /)
      
      real, parameter, dimension(8) :: Gqu_A_8 = &
      (/ 0.1012285362903763, &
         0.2223810344533745, &
         0.3137066458778873, &
         0.3626837833783620, &
         0.3626837833783620, &
         0.3137066458778873, &
         0.2223810344533745, &
         0.1012285362903763 /)
         
!     real, parameter, dimension(9) :: Gqu_t_9 = & 
!     (/ -0.9681602395076261, &
!        -0.8360311073266358, &
!        -0.6133714327005904, &
!        -0.3242534234038089, &
!         0.0000000000000000, &
!         0.3242534234038089, &
!         0.6133714327005904, &
!         0.8360311073266358, &
!         0.9681602395076261 /)
!     
!     real, parameter, dimension(9) :: Gqu_A_9 = & 
!     (/  0.0812743883615744, &
!         0.1806481606948574, &
!         0.2606106964029354, &
!         0.3123470770400029, &
!         0.3302393550012598, &
!         0.3123470770400029, &
!         0.2606106964029354, &
!         0.1806481606948574, &
!         0.0812743883615744 /)


!     real, parameter, dimension(14) :: Gqu_t_14 = & 
!     (/ -0.9862838086968123, &
!        -0.9284348836635735, &
!        -0.8272013150697650, &
!        -0.6872929048116855, &
!        -0.5152486363581541, &
!        -0.3191123689278897, &
!        -0.1080549487073437, &
!         0.1080549487073437, &
!         0.3191123689278897, &
!         0.5152486363581541, &
!         0.6872929048116855, &
!         0.8272013150697650, &
!         0.9284348836635735, &
!         0.9862838086968123 /)
!     
!     real, parameter, dimension(14) :: Gqu_A_14 = & 
!     (/  0.0351194603317519, &
!         0.0801580871597602, &
!         0.1215185706879032, &
!         0.1572031671581935, &
!         0.1855383974779378, &
!         0.2051984637212956, &
!         0.2152638534631578, &
!         0.2152638534631578, &
!         0.2051984637212956, &
!         0.1855383974779378, &
!         0.1572031671581935, &
!         0.1215185706879032, &
!         0.0801580871597602, &
!         0.0351194603317519 /)
         
      end module Gquadrature 
      
      module soilVars
!     use gloVars, only: dp
      save
      integer ::  N !number of layers. HALF-SPACE at N+1
      real*8, dimension(:),  allocatable :: Z,RHO,BETA0,ALFA0
      complex*16 ,dimension(:), allocatable :: ALFA,BETA,AMU,LAMBDA
      real*8 :: Qs,Qp
      end module soilVars
      
      module sourceVars
      integer :: efsource
      real    :: xfsource,zfsource
      real*8  :: nxfsource,nzfsource
      logical :: intfsource
      end module sourceVars
      
      module waveNumVars
      !frequency loop vars:
      integer,save      :: NFREC
      real*8   ,save    :: FREC,DFREC,OME,OMEI
      !Discrete Wave-number:
      real   ,save      :: DK    ! delta k on discrete wave number
      integer,save      :: NMAX
      real   ,save      :: LFREC,Qq
      complex*16, dimension(:), allocatable :: expK
      complex*16, dimension(:,:), allocatable :: Vout
      real   ,save      :: minKvalueW, minKvalueU
      integer,save      :: lapse,KtrimStart
      end module waveNumVars
      
      module refSolMatrixVars
      !reference solution variables:
      complex*16 :: cOME  
      complex*16, save, allocatable :: B(:,:)
      complex*16, dimension(:,:,:), allocatable :: Ak
      complex*16, dimension(:), allocatable :: this_B                             
      integer, dimension(:), allocatable :: IPIV
      integer :: info    
      end module refSolMatrixVars
      
      module waveVars
!     save
      real, save :: Escala
      real, save :: Theta !grados
      real*8, save :: Dt,maxtime  !segundos
      real, save :: t0
      integer, save :: tipoPulso ! 0 dirac; 1 ricker
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
      real*8,save, dimension(:,:,:), allocatable :: Nodes_xz_n
      real,save, dimension(:,:), allocatable :: Xcoord 
            
      ! length and normals at midpoint of BigSegments:
      real, allocatable :: midPoint(:,:)
      real,save, allocatable :: lengthXI(:),layerXI(:)
      real*8,save, allocatable :: normXI(:,:)
      end module GeometryVars
      
      module resultVars
       type FFres
         complex*16 :: W,U,Tx,Tz
       end type FFres
      
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
      !complex*16, dimension(:,:,:,:), allocatable :: qlmk
       logical :: shouldUseQuadrature
       end type tipo_qlmk
       
       type Punto
        !                 ,--- 1 for x, 2 for z
        real    :: center(2)
        real    :: bord_A(2),bord_B(2)
        real*8  :: normal(2)
        real    :: length
        integer :: layer
        logical :: isBoundary 
        logical :: isOnInterface
        logical :: guardarFK
        logical :: guardarMovieSiblings
        integer :: boundaryIndex
      !                       ,--- f: 1...nfrec+1
      !                       | ,--- k: 1...NMAX+1 / 2*NMAX
      !                       | | ,--- iMec: 1:5
      !                       | | | 
        complex*16, dimension(:,:,:), allocatable :: FK
        complex*16, dimension  (:,:), allocatable :: FKh,FKv
        
      ! para las funciones de Green:
      !                            ,--- xi
      !                            |  
!       type(tipo_qlmk), dimension(:), allocatable     :: GT_k
        
      !               (imec) l ---, ,--- m (1 = horz, 2= vert)
      !                  iGQ ---, | |   ,--- k
      !                 xi ---, | | |   | 
      !                       | | | |   | 
!       complex*16, dimension(:,:,:,:,  :), allocatable :: GT_gq_k
!       complex*16, dimension(:,:,:,:)    , allocatable :: GT_gq
        complex*16, dimension(:  ,:,:)    , allocatable :: GT_gq
!       complex*16, dimension(:,:,:,:,:,:,:), allocatable :: GT_gq_mov_k
!       complex*16, dimension(:,:,:,:,:)  , allocatable :: GT_gq_mov
        complex*16, dimension(:  ,:,:,:)  , allocatable :: GT_gq_mov
      !                               |
      !                               '--- iMov
        
        
        ! espectro campo total inquirePoints : 
      !                         ,--- f: 1...nfrec+1
      !                         | ,--- iMec: 1:2
        complex*16, dimension  (:,:), allocatable :: W 
        ! espectro campo total moviePoints : 
        complex*16, dimension(:,:,:), allocatable :: WmovieSiblings
      !                       |
      !                       `--- sibling index
      
      !               ,--- xXx (indice punto integracion Gaussiana)
      !               | ,--- (1,2) -> (x,z)
      !               | |
      real, dimension(:,:), allocatable :: Gq_xXx_coords
      real, dimension(:), allocatable :: Gq_xXx_C
      
       end type Punto   
      ! bondary elements:     ,--- POINT index / x (receptor)
      type (Punto), dimension(:), allocatable, save, target :: allpoints
      type (Punto), dimension(:), allocatable, save :: inqPoints
      type (Punto), dimension(:), allocatable, save :: moviePoints
      type (Punto), dimension(:), allocatable, save, target :: BouPoints !xi
      type (Punto), save, target :: Po
      integer, save :: nIpts, nMpts, NxXxMpts, nBpts, nPts,&
                       iPtini,iPtfin,mPtini,mPtfin,bPtini,bPtfin
      
      !                     ,-,--- nBpts x nBpts (considering 1 direction)
      complex*16, dimension(:,:), allocatable :: ibemMat
      complex*16, dimension(:), allocatable :: trac0vec 
      integer, dimension(:), allocatable :: IPIVbem
      
      !                                          ,--- ix
      !                                          | ,--- iz
      !                                          | | ,--- iMec
      !                                          | | | ,--- foto{2*Nfrec}
      complex*16, allocatable, save :: Sismogram(:,:,:,:)
      complex*16, dimension(:), allocatable :: auxKvect
      real*8,dimension(:), allocatable :: Hf,Hk !filtros 
      
      !type pointerTable
      !  real :: z
      !  integer :: e
      !  type (Punto) :: Pt(:)
      !  integer :: nXs
      !end type pointerTable
      !type (pointerTable), pointer :: PoTa(:),auxPoTa(:)
      integer, allocatable, dimension(:,:), save :: fixedPoTa,pota,auxpota
      integer :: nZs

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
!      real, dimension(4,2), save :: colorBounds
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
      
      module hank
      contains
      
      SUBROUTINE HANKELS(Z,H0,H1)
!     Z = COMPLEX ARGUMENT
!
!     COMPUTE SECOND KIND HANKEL FUNCTIONS H0 AND H1
!
      COMPLEX*16 :: Z,H0,H1,C,A,E,E2,ZH,P
      real*8 :: X,Y,R,PHI,J,AR
      X=REAL(Z)
      Y=AIMAG(Z)
      R=SQRT(X*X+Y*Y)
      PHI=ATAN2(Y,X)
      IF(R.LE.10.0)GO TO 20
      J=2.0*R
      C=(0.0,0.1250)/Z
      K=2
      P=C*C
      A=4.5*P
      P=7.5*P
      H0=1.0+C+A
      H1=1.0-3.0*C-P
10    I=4*K
      K=K+1
      DI=I
      DK=K
      A=A*C*(DI+1.0/DK)
      P=P*C*(DI-3.0/DK)
      H0=H0+A
      H1=H1-P
      AR=ABS(REAL(P))+ABS(AIMAG(P))
      IF(AR.GT.1.E-16.AND.K.LT.J)GO TO 10
      AR=0.785398163397448-X-PHI/2.0
      E=0.0
      IF(Y.GT.-160.0) E=0.7978845608028650/SQRT(R)*EXP(Y)*CMPLX(COS(AR),SIN(AR),8)
!     IF(X.EQ.0.0)E=CMPLX(0.0,AIMAG(E))
      IF(abs(X) .lt. 0.00001)E=CMPLX(0.0,AIMAG(E),8)
      H0=H0*E
      H1=H1*E*(0.0,1.0)
      GO TO 23
20    ZH=Z/2.0
      C=-ZH*ZH
      E=CMPLX(0.0,0.3183098861837910,8)
      E2=E*2.0
      A=1.0-E2*(0.5772156649015330+LOG(R/2.0))+PHI*0.636619772367582
      P=1.0
      K=1
      H0=A
      H1=A+E*(1.0-1.0/C)
25    A=A+E2/K
      P=P*C
      H0=H0+A*P
      K=K+1
      P=P/(K*K)
      H1=H1+(A*K+E)*P
      IF(ABS(REAL(P))+ABS(AIMAG(P)).GT.1.E-16)GO TO 25
      H1=H1*ZH
!     IF(X.NE.0.0)GO TO 23
      IF(abs(X) .gt. 0.00001)GO TO 23
      H0=CMPLX(0.0,AIMAG(H0),8)
      H1=CMPLX(REAL(H1),0.0,8)
23    RETURN
      END SUBROUTINE HANKELS
      end module hank
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
        do while( (abs(bf) > errortol) .and. (abs(a-b) > errorTol ) )
!       do while((.not.(bf == real(0,8) )) .and. (abs(a-b) > errorTol ))
!          print*,"go_",i
          if ((abs(af-cf) > errorTol) .and. (abs(bf-cf)>errorTol)) then
!          if ((af /= cf) .and. (bf /= cf)) then
          ! interpolación cuadrática inversa
            s = a * bf * cf / (af-bf) / (af-cf) + & 
            b*af*cf/(bf-af)/(bf-cf)+cx*af*bf/(cf-af)/(cf-bf)
          else
          ! regla de la secante
            s = b - bf * (b-a)/(bf-af)
          end if
          tmp2 = (3.0*a + b)/4.0
          if ( (.not.(((s > tmp2) .and. (s < b)) .or. & 
          ((s < tmp2) .and. (s > b))) ) .or. &
          (mflag .and. ((abs(s-b)) .ge. (abs(b-cx)/2.0 ))) .or. &
          ((.not. (mflag)) .and. ((abs(s-b)) .ge. (abs(cx-d)/2.0 ))) ) then
            s = (a+b) / 2.0
            mflag = .true.
          else
            if ((mflag .and. (abs(b-cx)< errorTol)) .or. &
            ((.not. (mflag)) .and. (abs(cx-d) < errorTol))) then
              s = (a+b) / 2.0
              mflag = .true.
            else
              mflag = .false.
            end if
          end if
           sf = polyVal(surf_poly0,degree,s)
           d = cx
           cx = b
           cf = bf
!          if (af * sf < real(0,8)) then
           if (af * sf < errorTol) then
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
             b = real(0.123456789)
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
!     integer, parameter             :: dp = selected_real_kind(15, 307)
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
      real*8, allocatable :: normalto(:,:) !function result value
      real*8, parameter      :: a = real(1,8) ! normal poing up parameter
!     real*8, parameter      :: tolerance = real(0.001,8) !to tell if its vertical
      real*8, dimension(degree) :: fprime !surf_poly derivative coeficients
      real*8 :: fprimeX, x0, mag 
      
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
      complex,    dimension(:), allocatable :: y
      real,       dimension(:), allocatable :: x
      real maxY,minY,maxYc,minYc,xstep,ystep
      integer :: i
      character(LEN=100) :: dumb
      CHARACTER(LEN=6)  :: CBUF
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      allocate(x(n))
      allocate(y(n))
      DO i = 1,n
        x(i) = Df*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
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
      CALL axspos (int(370,4) ,int(H+100,4)) !Axis system pos. Lower left corner
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
      
!     print*,'plotted ',trim(titleN)
!     !print*,''
      deallocate(x,y)
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
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      
!     allocate(x(n))
!     print*,size(y)
      DO i = 1,n
        x(i) = Dt*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
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
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
!     
      
!     allocate(x(n))
!     print*,size(y)
      DO i = 1,n
        x(i) = Dt*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
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
      complex*16, dimension(NFREC+1,NMAX,2),intent(in) :: thisFK
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
      
      
      if(verbose>=2)write(outpf,'(a,a,a)') "will plot ",trim(tt),"..."
            
      M = 0
      nombre(1)= '_w.png'
      nombre(2)= '_u.png'
      nombre(3)= '_s33.png'
      nombre(4)= '_s31.png'
      nombre(5)= '_s11.png'
      
      do i=1,NFREC
         vHorz(i) = 0 + (i-1) * real(DFREC,4) * 2 * pi ! frec angular
      end do
      !
      do ik=1,nmax
         vVert(ik) = 0 + (ik-1) * real(DK,4)
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
      do iMec = 1,2 !min(2,imecMax) !podrian ser más pero meh      
      write(titulo,'(a,a,I0,a,I0,a,I0,a,I0,a,a)') trim(tt),'[', &
      int(coords(1)),'.',abs(int((coords(1)-int(coords(1)))*10)),';', & 
      int(coords(2)),'.',abs(int((coords(2)-int(coords(2)))*10)),']', & 
      trim(nombre(iMec))
      
      if (form .eq. 1) then
       M = real(thisFK(1:NFREC,1:NMAX,iMec),4)
      elseif (form .eq. 2) then
       M = real(aimag(thisFK(1:NFREC,1:NMAX,iMec)),4)
      elseif (form .eq. 3) then
       M = real(log(1. + exp(p)*abs(thisFK(1:NFREC,1:NMAX,iMec))) / & 
           log(exp(p)+1.),4)
       M = M / maxval(M)
      end if
!     print *,titulo
      miV = minval(M)!;print *, "min",miV
      maV = maxval(M)!;print *, "max",maV
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
        maxY = max(maxval(Xcoord(:,2))+maxval(Xcoord(:,2))/10,real(Z(N+1),4))
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
      use soilVars!, only : N
      use debugStuff
      use resultVars
      use GeometryVars!, only : numberOffixedBoundaryPoints
      use ploteo10pesos
      use sourceVars
      use wavelets
      use MeshVars
!     use Gquadrature, only: Gquad_n
!     use waveVars, only : Uo!,Escala
      
      implicit none
      
      ! output direction   6 : screen      101: log file
      integer, parameter :: PrintNum = 6
      !Loop Counters
      integer :: J  ! frequency loop
      integer :: l,m,iP,iPxi,iP_x,i,ik,iz !small loop counter
      integer :: status 
      CHARACTER(len=400) :: path
      character(LEN=100) :: titleN
      character(LEN=10) :: tt
      character(LEN=9)  :: xAx,yAx
      real*8 :: factor,k
      character(10) :: time
      real :: startT,finishT,tstart,tend
      logical :: thereisavirtualsourceat
      call cpu_time(startT)
      call cpu_time(tstart)
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
      call system('mkdir outs')
      CALL chdir("outs",status)
      if (status .eq. 0) call system("rm *.*")
      if (status .eq. 0) call chdir("..")
      
      call getMainInput
      
      call getSoilProps (PrintNum)    
      nIpts=0; nMpts=0; nBpts = 0
      iPtfin = 0
      mPtfin = 0
      KtrimStart = 1
      allocate(vout(1,2))
      if (getInquirePointsSol) call getInquirePoints(PrintNum)
      call getVideoPoints(PrintNum)
      if (workBoundary) call getTopography(PrintNum)
      NPts = nIpts + nMpts
      write(PrintNum,'(a,I0)') 'Number of fixed receptors: ',npts
      allocate (allpoints(Npts))
      allocate (Ak(4*N+2,4*N+2,nmax+1))
      ALLOCATE (B (4*N+2,nmax+1)) ! una sola fuente
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
          allocate(allpoints(iP)%FK(NFREC+1,NMAX,2))
          allpoints(iP)%FK = 0
      end if; end do!
      
      do iP=iPtini,iPtfin !solo inquirepoints:
         allocate(allpoints(iP)%W(NFREC+1,imecMax)); allpoints(iP)%W = 0
      end do
      
      if (makeVideo) then; do iP=mPtini,mPtfin
         allocate(allpoints(iP)%WmovieSiblings(NxXxMpts,NFREC+1,imecMax))
         allpoints(iP)%WmovieSiblings = 0
      end do; end if!
      
      if (workBoundary) then
      ! G para resolver el campo difractado por topografía
       do iP_x = iPtini,iPtfin
         allocate(allpoints(iP_x)%GT_gq(nBpts,2,2))
       end do
       if (makeVideo) then
         do iP_x = mPtini,mPtfin
           allocate(allpoints(iP_x)%GT_gq_mov(nBpts,2,2,NxXxMpts))
         end do
       end if
      end if
      
      ! source application point:
      call getsource
      
      if (verbose .ge. 1) then
        write(PrintNum,'(a,F8.2,a,F8.2,a)') & 
                            "source is at (",xfsource,",",zfsource,")"
      end if
      
      call makeTaperFuncs_cutF_cutK(20,0.8,0.8)
      call expPIK(expK)
      call waveAmplitude(PrintNum)
      
      call preparePointerTable(.true.,PrintNum)
      call cpu_time(tend)
      print*,"Pre-prosses took ",tend-tstart,"seconds"
      
      write(PrintNum,'(A)', ADVANCE = "NO") & 
                                 "J |             omega            |    [Hz]   "
      if(WORKBOUNDARY)then
      write(PrintNum,'(A)', ADVANCE = "YES") " | minS WL [m] | Nel | t span [s]"
      else 
          write(PrintNum,'(A)', ADVANCE = "YES") " | t span"
      end if
      DO J=1,NFREC+1
      ! complex frecuency for avoiding poles and provide damping
        FREC=DFREC*real(J-1) ! Hz
        if (J .eq. 1)  FREC = DFREC*0.01
        call cpu_time(tstart)
        OME=2.0*PI*FREC !rad/s
        !periodic sources damping
        cOME = cmplx(OME, - DFREC / periodicdamper,8) !* cmplx(1.0, -1.0/2.0/Qq,8)
         
        ! Azimi attenuation
        do l=1,N+1
          Alfa(l) = Alfa0(l)*(cmplx(1+1/(pi*Qp)*log(ome/(2*pi)),1/(2*Qp),8))
          beta(l) = beta0(l)*(cmplx(1+1/(pi*Qs)*log(ome/(2*pi)),1/(2*Qs),8))
          aMU(l) = RHO(l) * beta(l)**2
          Lambda(l) = RHO(l)*Alfa(l)**2 - real(2.)*aMU(l)
        end do
         
!       ! complex frecuency for avoiding poles and provide damping
!       FREC=DFREC*real(J-1)
!       if (J .eq. 1) FREC = DFREC*0.001 
!       OME=2.0*PI*FREC
!       OMEI=0.7*PI/TW
!       OMEI=DFREC / 2.0
        ! Bouchon (2003) OMEI entre -pi/T y -2pi/T ; T= 2pi/DFREC
!       COME=CMPLX(OME, -OMEI*1.0) !periodic sources damping
!       COME=COME*(UR - UI/2.0/Qq) !histeretic damping
!       cOME = cmplx(OME, -DFREC / periodicdamper,8) * cmplx(1.0, -1.0/2.0/Qq,8)
        
        if(verbose .eq. 1)then
 write(PrintNum,'(I0,A,EN13.2,1x,EN13.2,A,EN10.2,a)', ADVANCE = "NO") &
 J,' | ',REAL(COME),aimag(COME),'i | ',FREC,' | '
        else if (verbose .gt. 1) then
 write(PrintNum,'(A,I0,A,EN13.2,1x,EN13.2,A,EN18.2)', ADVANCE = "YES") &
 'w(',J,') ',REAL(COME),aimag(COME),'i | ',FREC
        end if 
        
      ! Subsegment the topography if neccesssary
      if (workBoundary) then 
         call subdivideTopo(J,PrintNum)
!        call allocintegPoints(J)
         call preparePointerTable(.false.,PrintNum)
      end if
         
      Ak = 0
      Do ik=1,nmax+1 !WAVE NUMBER LOOP-
         k = real(ik-1) * dK
         if (ik .eq. 1) k = dk * 0.001
         call matrixA_borderCond(Ak(:,:,ik),k,cOME,PrintNum)
         call inverseA(Ak(:,:,ik),4*N+2)
      end do ! ik
      
      if (workboundary) then
        trac0vec = 0
        ibemMat = 0
      end if
        ! for the discrete wave number
!     LFREC=2000 + L*exp(-PI*(FLOAT(J-1)/NFREC)**2)  !EMPIRICAL (improve)
!     DK=2.0*PI/LFREC
      
!     NMAX=(OME/minval(BETA))/DK+1                   !EMPIRICAL (improve)
!     NMAX=2*NMAX+100                                !EMPIRICAL (improve)
      
      call crunch(0,J,cOME,PrintNum)
      !           ^--- the real source
      if (workboundary) then
      if (verbose .ge. 2) then
        call showMNmatrixZabs(2*nBpts,1, trac0vec ,"tindp",PrintNum)
      end if
!     stop !stop "only the vector"
      do iz = 1,nZs
        if (thereisavirtualsourceat(iz)) then
          call crunch(iz,J,cOME,PrintNum)
        end if 
      end do !iz
      
      if (verbose .ge. 2) then
        call showMNmatrixZ(2*nBpts,2*nBpts,ibemMat,"ibMat",PrintNum)
      end if
!     stop "nice"
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
         CALL chdir("..")
      end if! verbose 1
      
      if(verbose .ge. 2) write(PrintNum,'(a)') "add diffracted field by topography"
      do iP_x = iPtini,iPtfin !cada receptor X
        do iPxi = 1,2*nBpts,2 !cada fuente virtual (dos direcciones)
          do l=0,1 !componete desp. en receptor X  (x,z)
            do m=0,1 !direc. fuente
        allpoints(iP_x)%W(J,2-l) = allpoints(iP_x)%W(J,2-l) + &
        allpoints(iP_x)%GT_gq(ceiling(iPxi/2.),2-l,m+1) * trac0vec(iPxi+m)
!       integralEq16(iP_x,l,ceiling(iPxi/2.),m) * trac0vec(iPxi+m)
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
        allpoints(iP_x)%WmovieSiblings(iP,J,2-l) = & 
        allpoints(iP_x)%WmovieSiblings(iP,J,2-l) + &
        allpoints(iP_x)%GT_gq_mov(ceiling(iPxi/2.),2-l,m+1,iP) * trac0vec(iPxi+m)
!       integralEq16mov(iP,iP_x,l,ceiling(iPxi/2.),m) * trac0vec(iPxi+m)
            end do !m
          end do !iPxi
        end do ! l
       end do ! iP
      end do !iP_x
      end if
            
      end if !workboundary
      call cpu_time(tend)
      if (verbose .eq. 1) write(PrintNum,'(E9.1)') tend-tstart
      if (verbose .gt. 1) write(PrintNum,'(a,E9.1)') "time span = ",tend-tstart
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
           call system("mkdir outs")
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
      use wavenumvars, only : minKvalueW, minKvalueU,lapse
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
      read(35,*) periodicdamper; print*,"Periodic sources damping factor ", periodicdamper
      read(35,*)
      read(35,*) minKvalueW; print*,"min k para despreciar W= ", minKvalueW
      read(35,*) minKvalueU; print*,"min k para despreciar W= ", minKvalueU
      read(35,*) lapse; print*,"number of minimum values commited ", lapse
      read(35,*)
      read(35,*) plotStress
      read(35,*)
      read(35,*) plotModule
      close(35)
      if (plotStress) then 
       imecMax = 4 ! dos desps , dos tracs
      else 
       imecMax = 2 ! dos desps
      end if
      plusOne = 0
      if (plotModule) plusOne = 1
      end subroutine getMainInput
      subroutine getSoilProps (outpf)
      use soilVars
      use waveNumVars
      use waveVars, only : dt,t0
      use gloVars, only : verbose!,PI
      implicit none
      integer, intent(in) :: outpf
      logical :: lexist
      real :: H, ALF, BET, RO!, rKmax
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
      ALLOCATE (BETA0(N+1))
      ALLOCATE (ALFA0(N+1))
      ALLOCATE (LAMBDA(N+1))
      ALLOCATE (RHO(N+1))
      
      Z(1)=real(0,8)
      if (verbose >= 1) then
       write(outpf,'(a)')& 
              '        depth       alpha0    beta0      mu0     rho      lambda0'
      end if
      DO J=1,N
         READ(7,*) H, ALF, BET, RO
         Z(J+1)=Z(J)+real(H)
         AMU(J)=cmplx(RO*BET**2,0.0,8)
         BETA0(J)=BET
         ALFA0(J)=ALF
         RHO(J) = RO
         LAMBDA(J)=cmplx(RHO(J)*ALF**2 - real(2.)*real(AMU(J)),0.0,8)
!         BEALF=SQRT((0.5-ANU)/(1.0-ANU)) !IF POISSON RATIO IS GIVEN
!         ALFA(J)=BET/BEALF
       if (verbose >= 1) then
         write(outpf,& 
       '(F7.1,A,F7.1,2x, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2)') & 
       Z(J),' - ',Z(J+1),ALFA0(J),BETA0(J),real(AMU(J)),  RHO(J), real(LAMBDA(J))
       end if
      END DO
      
      READ(7,*) H, ALF, BET, RO 
      AMU(N+1)=cmplx(RO*BET**2,0.0,8)
      BETA0(N+1)=BET
      ALFA0(N+1)=ALF
      RHO(N+1) = RO
      LAMBDA(N+1)=cmplx(RHO(n+1)*ALF**2 - real(2.)*real(AMU(n+1)),0.0,8)
      if (verbose >= 1) then
!      write(outpf,'(F7.1,A,2x,F8.2,3x,F7.1,6x,F7.1)') &
!      Z(N+1),' -    inf.',AMU(N+1),BETA0(N+1),ALFA0(N+1)
       write(outpf,'(F7.1,A,2x,F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2)') & 
       Z(N+1),' -   inf. ',ALFA0(N+1),BETA0(N+1),real(AMU(N+1)),RHO(N+1),real(LAMBDA(N+1))
       write(outpf,'(a)')' '
      end if
!      BEALF=SQRT((0.5-ANU)/(1.0-ANU))   !IF POISSON RATIO IS GIVEN
!      ALFA(N+1)=BET/BEALF
      
      READ(7,*)
      READ(7,*)DFREC,NFREC,DK,nmax,Qq,t0 
      close(7)
      Qs = Qq
      Qp = Qq
      ! aplicamos el amortiguamiento en las velocidades
      do J=1,N+1
        !non dispersive           
        Alfa(J) = Alfa0(J) * (cmplx(1.0,1/(2*Qq)))
        beta(J) = beta0(J) * (cmplx(1.0,1/(2*Qq)))
        aMU(J) = RHO(J) * beta(J)**2
        Lambda(J) = RHO(J)*Alfa(J)**2 - real(2.)*aMU(J)
        ! or Azimi every frec -> OK
      end do
      
      ! decidimos el máximo número de onda a partir
      ! de la frecuencia máxima y la velocidad de propagación
      ! mínima en el medio mas 20% por si a caso.
!     rKmax = (DFREC * NFREC)/(minval(beta))
!     rKmax = rKmax/DK
!     nMAX = ceiling(rKmax)
      
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
!     use soilVars, only : Z,N
!     use gloVars, only : dp
      implicit none
      integer, intent(in) :: outpf
      integer :: i,thelayeris
      logical :: lexist, tellisoninterface
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
      
      !encontrar el layer en el que estan o 0 si está sobre la interfaz
!          print*,"z=",nonBoundPoints(1)%center(i,2),"....."
           inqPoints(i)%layer = thelayeris(real(inqPoints(i)%center(2)))
           
           if (auxGuardarFK .eq. 1 ) then
              inqPoints(i)%guardarFK = .true.
           end if
           
           inqPoints(i)%isOnInterface = tellisoninterface(real(inqPoints(i)%center(2)))
      end do
!     stop 0
      end subroutine getInquirePoints   
      
      function thelayeris(zi)
      use soilVars, only : Z,N
      implicit none
      integer :: e,thelayeris
      real :: zi
      do e=1,N
         if(real(zi) .lt. Z(e+1)) then
            exit
         end if
      end do
      thelayeris = e
      end function thelayeris
      
      function tellisoninterface(zi)
      use soilVars, only : Z,N
      implicit none
      real ::  errT = 0.001
      integer :: e
      real :: zi
      logical :: tellisoninterface
      tellisoninterface = .false.
      do e=1,N+1
        if((Z(e)-errT .lt. real(zi)) & 
             .and. (real(zi) .lt. Z(e)+errT)) then
           tellisoninterface = .true.
        end if
      end do
      
      end function tellisoninterface
      subroutine getsource
      use wavevars, only: Escala,Ts,Tp, tipoPulso,maxtime
      use sourceVars
      use resultvars, only:Po
      implicit none
      integer :: thelayeris
      logical :: lexist, tellisoninterface
      
      inquire(file="source.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="source.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "source.txt" on Working directory' 
      end if
      READ(77,*)      
      READ(77,*) xfsource, zfsource, nxfsource, nzfsource
      READ(77,*)      
      READ(77,*) Escala
      READ(77,*)
      READ(77,*) tipoPulso
      READ(77,*) 
      READ(77,*) Ts
      READ(77,*)
      READ(77,*) Tp
      READ(77,*)
      READ(77,*) maxtime
      close(77)
      
      efsource = thelayeris(real(zfsource))
      intfsource = tellisoninterface(real(zfsource))
      
      Po%center(1) = xfsource
      Po%center(2) = zfsource
      Po%normal(1) = nxfsource
      Po%normal(2) = nzfsource
      Po%layer = efsource
      Po%isOnInterface = intfsource
      allocate(Po%Gq_xXx_coords(1,2))
      Po%Gq_xXx_coords(1,1) = xfsource
      Po%Gq_xXx_coords(1,2) = zfsource
      end subroutine getsource
      subroutine getVideoPoints(outpf)
      use resultVars, only : moviePoints, nMpts, & 
                             iPtfin,mPtini,mPtfin
      use meshVars
      use gloVars,only:verbose,PI,makevideo
      use soilVars, only : Z,N
      use waveNumVars, only : DK,NMAX
      use resultvars, only: NxXxMpts
!     use wavevars, only: Escala,Ts,Tp, tipoPulso,maxtime
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
!     READ(7,*)
!     READ(7,*)
!     READ(7,*) Escala
!     READ(7,*)
!     READ(7,*) tipoPulso
!     READ(7,*) 
!     READ(7,*) Ts
!     READ(7,*)
!     READ(7,*) Tp
!     READ(7,*)
!     READ(7,*) maxtime
      close(7)
      
      if (makeVideo .eqv. .false.) return 
        
        write(outpf,'(a,F12.2,a)') "Lfrec = 2*pi/DK = ",2*PI/DK, "m"
        DeltaX = pi / (nMax * DK)
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
!       if(boundingZn == 0.0) then
        if(abs(boundingZn)<errT) then
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
                         MeshDXmultiplo,"] = ", NxXxMpts," puntos"
          write(outpf,'(a,F9.1,a,F9.1,a,F10.3,a,I0,a)') & 
                     "z :[",boundingZn," ... ",boundingZp, & 
                     "] :: ",MeshDZ,"m      = ",nz," puntos"
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
            moviePoints(iP)%normal(1) = 0.0
            moviePoints(iP)%normal(2) = 1.0
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
      use soilVars, only : Z,N!,BETA
      use waveNumVars, only : NMAX,vout
      use ploteo10pesos
      use resultVars, only : BouPoints, nBpts, & 
                             mPtfin,bPtini,bPtfin,iPtfin, &
                             ibemMat,trac0vec,IPIVbem 
!     use refSolMatrixVars, only: Bb
      use Gquadrature, only: Gqu_n => Gquad_n, & 
                             Gqu_t => Gqu_t_3, & 
                             Gqu_A => Gqu_A_3
      
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
      real*8 :: bigN
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
!         if (anint(y(iXI) * 1000) == anint(Z(e) * 1000)) then
          if (abs(anint(y(iXI) * 1000) - anint(Z(e) * 1000)) < errT) then
           if (verbose >= 2) then; write(outpf,*)"equal"; end if
          else 
      ! si es un segmento que se forma recorriendo los puntos hacia Z+
            if (y(iXI) < Z(e)  .AND. y(iXI+1) > Z(e)) then
               !hay un cambio de estrato
               nuevoPx = splitatY(surf_poly,degree,real(Z(e),4),x(iXI),x(iXI+1))
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
               nuevoPx = splitatY(surf_poly,degree,real(Z(e),4),x(iXI),x(iXI+1))
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
      ! reuglarizar normales
      bigN = 100000
      do iXi = 1,nXi-1
        normXI(iXi,1) = anint(normXI(iXi,1)* bigN)/bigN
        normXI(iXi,2) = anint(normXI(iXi,2)* bigN)/bigN
      end do
      
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
      BouPoints(:)%layer = int(layerXI)
      BouPoints(:)%isBoundary = .true.
      BouPoints(:)%isOnInterface = isOnIF
      BouPoints(:)%guardarFK = .false.
      BouPoints(:)%guardarMovieSiblings = .false.
      
      do iX=1,nBpts !van vacías porque esto cuenta para cada frecuencia
!       allocate(BouPoints(iX)%FKh(NMAX+1,5)); BouPoints(iX)%FKh = 0
!       allocate(BouPoints(iX)%FKv(NMAX+1,5)); BouPoints(iX)%FKv = 0
        allocate(BouPoints(iX)%FK(1,2*NMAX,5)) 
        allocate(BouPoints(iX)%W(1,2))
        
      ! Coordenadas/pesos de integración Gaussiana:
        allocate(BouPoints(iX)%Gq_xXx_coords(Gqu_n,2))
        allocate(BouPoints(iX)%Gq_xXx_C(Gqu_n))
        
        BouPoints(iX)%boundaryIndex = iX
        
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
!       allocate(BouPoints(iX)%GT_k(nBpts))
        
!       allocate(BouPoints(iX)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); BouPoints(iX)%GT_gq_k = 0
!       allocate(BouPoints(iX)%GT_gq(nBpts,Gqu_n,5,2)); BouPoints(iX)%GT_gq = 0
!       allocate(BouPoints(iX)%GT_gq(nBpts      ,5,2)); BouPoints(iX)%GT_gq = 0
        
      end do !iX
      
      
      deallocate(Vout)
      allocate(Vout(2*nBpts,2))
      allocate(ibemMat(2*nBpts,2*nBpts))
      allocate(trac0vec(2*nBpts))
      allocate(IPIVbem(2*nBpts))
      
      end subroutine getTopography
      subroutine preparePointerTable(firstTime,outpf)
      use resultVars, only : nPts,allpoints,nBpts,BouPoints,auxpota,pota,fixedpota,nZs,Punto
      use debugstuff
      use Gquadrature, only : Gquad_n
      use glovars, only : verbose,workBoundary
      ! tabla con las coordenadas Z (sin repetir).
      implicit none
      logical,intent(in) :: firstTime
      integer,intent(in) :: outpf
      integer :: i,j,iP
      logical :: nearby,esnueva
      type(Punto), pointer :: PX
      ! si es la primera vez que corre sólo agregamos los allpoints 
      if (firstTime) then
       if (verbose .ge. 2) write(outpf,*) "generating master table"
       allocate(auxpota(npts,npts+2))
       auxpota = 0
       ! siempre agregamos el primer punto [ allpoints(1)% ]
       nZs = 1
       auxpota(1,1) = 1 !nXs allpoints
!      auxpota(1,2) = 0 !nXs boupoints
       auxpota(1,3) = 1 !4,5,6 ... allpoints -7,-8,-9 ... boupoints
       !         '--- el (3) siempre está
       ! agregamos los demás sin repetir
       if (nPts>=2) then 
       do iP = 2,nPts
         esnueva = .true.
         do i = 1,nZs
           if (nearby(allpoints(iP)%center(2),allpoints(auxpota(i,3))%center(2),0.001)) then
             !agregar a coordenada Z existente
             auxpota(i,1) = 1 + auxpota(i,1)
             auxpota(i,auxpota(i,1) + 2) = iP
             esnueva = .false.
             exit
           end if
         end do
         if (esnueva) then
           !nueva coordenada Z
           nZs = nZs + 1
           auxpota(nZs,1) = 1
           auxpota(nZs,3) = iP
         end if
       end do
       end if
       j = maxval(auxpota(:,1))
       
       allocate(PoTa(nZs,2+j))
       Pota = auxpota(1:nZs,1:2+j)
       deallocate(auxpota)
       if (workBoundary) then
         allocate(fixedPoTa(nZs,2+j))
         fixedPota = Pota
       end if
      else !---------------------------------------------------------------------------------------
       ! Dada la tabla de los puntos fijos.
       ! Agregar los puntos gaussianos de los segmentos cercanos a la
       ! fuerza virtual [zf] y los centros de los segmentos restantes. 
       if (verbose .ge. 2) write(outpf,*) "updating master table"
       if (allocated(pota)) deallocate(pota)
       if (allocated(auxpota)) deallocate(auxpota)
!      print*,"max size:  ",nZs + nBpts * Gquad_n," x ",2 + maxval(fixedPota(:,1)) + nBpts * Gquad_n
       allocate(auxpota(nZs + nBpts * Gquad_n, 2 + maxval(fixedPota(:,1)) + nBpts * Gquad_n))
       auxpota = 0
!      j = maxval(fixedPota(:,1))
       nZs = size(fixedPota,1)
!      print*,"nzs=",nzs
!      print*,"fixed"
!      print*,size(fixedPota,1),size(fixedPota,2)
!      call showMNmatrixI(size(fixedPota,1),size(fixedPota,2),fixedPota,"fixed",outpf)
!      print*,""
!      print*,"auxpota"
!      print*,size(auxpota,1),size(auxpota,2)
       auxpota(1:size(fixedPota,1),1:size(fixedPota,2)) = fixedPota
!      read(5,*)
!      print*,"nBpts=",nBpts
       do iP = 1,nBpts
!        print*,"ip:",ip
         esnueva = .true.
         do i = 1,nZs
          ! diferenciar de que grupo es la coordenda
           if (associated(PX)) then 
              nullify(PX)
           end if!
           if (auxpota(i,3) .gt. 0) then
             PX => allpoints(auxpota(i,3))!; print*,"allp "
           else
             PX => boupoints(abs(auxpota(i,3)))!; print*,"boup "
           end if
          ! revisar si existe
           if (nearby(boupoints(iP)%center(2),PX%center(2),0.001)) then
!            print*,"estan cerca"
             if (auxpota(i,3) .lt. 0) then 
               if (.not.(nearby(PX%length,boupoints(iP)%length,0.01)) .or. &
              (.not.(nearby(real(abs(PX%normal(1)),4),real(abs(boupoints(iP)%normal(1)),4),0.01)))) then
                 esnueva = .true. ;print*," they are different"
                 exit
               end if
             end if
             !inscribir a coordenada Z existente
             auxpota(i,2) = 1 + auxpota(i,2)
             auxpota(i,auxpota(i,1) + auxpota(i,2) + 2) = - iP
             esnueva = .false.
!            print*,"inscrito. renglon:",i,"(",auxpota(i,1)," ",auxpota(i,2),")"
             exit
           end if
         end do ! i
         if (esnueva) then
           !nueva coordenada Z
           nZs = nZs + 1
           auxpota(nZs,1) = 0
           auxpota(nZs,2) = 1
           auxpota(nZs,3) = - iP
!          print*,"nuevo Z ahora son ",nZs
         end if
       end do ! iP
       j = maxval(auxpota(:,1) + auxpota(:,2))
       allocate(PoTa(nZs,2+j))
       Pota = auxpota(1:nZs,1:2+j)
       deallocate(auxpota)
      end if
      ! done
!     print*,"nZs=",nZs
      if (verbose .ge. 2) call showMNmatrixI(nZs,2+j,pota,"po_ta",outpf)
!     if (firstTime .eqv. .false.) stop "preparePointerTable"
      end subroutine preparePointerTable
      
      function nearby(a,b,bola)
      implicit none
      real, intent(in) :: a,b,bola
      logical :: nearby
      nearby = .false.
      if ((b-bola .le. a) .and. (a .le. b+bola)) then 
        nearby = .true.
      end if
!     print*, b-bola, "<=", a, "<=",b+bola ," :: ", nearby
      end function nearby 
      subroutine subdivideTopo(iJ,outpf)
      !Read coordinates of collocation points and fix if there are
      !intersections with inferfaces. Also find normal vectors of segments.
      use GeometryVars! x,y,Xcoord,surf_poly,nXI
      use gloVars!, only : verbose, makeVideo, getInquirePointsSol
      use fitting
      use soilVars, only : Z,N,BETA
      use waveNumVars, only : DFREC,NMAX,vout
      use resultVars, only : BouPoints, nBpts, & 
                             bPtini,bPtfin,iPtini,iPtfin, &
                             allpoints,nPts,mPtini,mPtfin,NxXxMpts, &
                             ibemMat,trac0vec,IPIVbem
!     use refSolMatrixVars, only: Bb
      use ploteo10pesos
      use Gquadrature, only: Gqu_n => Gquad_n, & 
                             Gqu_t => Gqu_t_3, & 
                             Gqu_A => Gqu_A_3

      implicit none
      integer, intent(in) :: outpf,iJ
      real :: deltaX,deltaZ,maxLen,leniXI
      integer :: J,iXI,iX,e,indice, AllocateStatus,idx
      character(LEN=100) :: txt
      real     :: errT = 0.001
      logical, dimension(:), allocatable  :: isOnIF
      real :: thisFrec
      real, dimension(2) :: midiXI
      logical :: smallEnough, huboCambios
      real, dimension(:), allocatable :: auxA,auxB
      
      real, dimension(2) :: norm_comp
      real, dimension(2) :: ABp !x or y coords of points A and B
      integer :: i,xory,xoryOtro
      real :: xA,yA,xB,yB
      real :: interpol, smallestBeta
      integer :: iNstart,iNfinish,iN
      
      huboCambios = .false.
      allocate(isOnIF(1))
      J = iJ
      if (J .lt. 2) J=2
      thisFrec =  real(DFREC,4)*real(J-1) ! Hz
      allocate(auxA(1))
      allocate(auxB(1))
      smallestBeta = minval(real(abs(beta(:)),4)) !smallest overall
      !------- subdivide to SEE the smallest wavelenght at every stratum 
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
           maxLen = real(abs(BETA(e)),4) / (multSubdiv * thisFrec)
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
      
      if (verbose .ge. 1) then
       CALL chdir("outs")
       write(txt,'(a,I0,a)') 'Division_at[J=',iJ,'].pdf'
       call drawGEOM(txt,.false.,outpf)
       CALL chdir("..")
       if (verbose .gt. 1) write(outpf,'(a,I0,a)')"Boundary was segmented, now it is made of: (",nXI-1,") elements"
       
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
        if(allocated(BouPoints(iX)%GT_gq)) deallocate(BouPoints(iX)%GT_gq)
!       do iXi=1,size(BouPoints(iX)%GT_k)
!         if(allocated(BouPoints(iX)%GT_k(ixi)%qlmk)) deallocate(BouPoints(iX)%GT_k(iXi)%qlmk)
!       end do
!       deallocate(BouPoints(iX)%GT_k)
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
      BouPoints(:)%layer = int(layerXI)
      BouPoints(:)%isBoundary = .true.
      BouPoints(:)%isOnInterface = isOnIF
      BouPoints(:)%guardarFK = .false.
      BouPoints(:)%guardarMovieSiblings = .false.
      
      ! para resolver las densidades de fuerza IBEM:
      do iX=1,nBpts !van vacías porque esto cuenta para cada frecuencia
        ! start saving space for the green functions:
!       allocate(BouPoints(iX)%GT_k(nBpts))
!       allocate(BouPoints(iX)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); BouPoints(iX)%GT_gq_k = 0
!       allocate(BouPoints(iX)%GT_gq(nBpts,Gqu_n,5,2)); BouPoints(iX)%GT_gq = 0
        allocate(BouPoints(iX)%GT_gq(nBpts      ,2,2)); BouPoints(iX)%GT_gq = 0
          
!       allocate(BouPoints(iX)%FKh(NMAX+1,5)); BouPoints(iX)%FKh = 0
!       allocate(BouPoints(iX)%FKv(NMAX+1,5)); BouPoints(iX)%FKv = 0
        allocate(BouPoints(iX)%FK(1,2*NMAX,5))  
        allocate(BouPoints(iX)%W(1,2))
        
        allocate(BouPoints(iX)%Gq_xXx_coords(Gqu_n,2))
        allocate(BouPoints(iX)%Gq_xXx_C(Gqu_n))
       
        BouPoints(iX)%boundaryIndex = iX
       
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
!       if(allocated(allpoints(ix)%GT_k)) deallocate(allpoints(ix)%GT_k) 
!       allocate(allpoints(ix)%GT_k(nBpts))
      
!       deallocate(allpoints(ix)%GT_gq_k)
!         allocate(allpoints(ix)%GT_gq_k(nBpts,Gqu_n,5,2,nmax+1)); allpoints(ix)%GT_gq_k = 0
      end do!
      do iX=iPtini,iPtfin 
          if(allocated(allpoints(ix)%GT_gq)) deallocate(allpoints(ix)%GT_gq)
!         allocate(allpoints(ix)%GT_gq(nBpts,Gqu_n,5,2)); allpoints(ix)%GT_gq = 0
          allocate(allpoints(ix)%GT_gq(nBpts      ,2,2)); allpoints(ix)%GT_gq = 0
      end do
      
      if (makeVideo) then
        do iX = mPtini,mPtfin
          if(allocated(allpoints(ix)%GT_gq_mov)) deallocate(allpoints(ix)%GT_gq_mov)
!         allocate(allpoints(ix)%GT_gq_mov(nBpts,Gqu_n,5,2,NxXxMpts)); allpoints(ix)%GT_gq_mov = 0
          allocate(allpoints(ix)%GT_gq_mov(nBpts      ,2,2,NxXxMpts)); allpoints(ix)%GT_gq_mov = 0
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
      if (verbose .eq. 1) write(outpf,'(EN11.1,a,I0,a)', ADVANCE = "NO") smallestBeta/thisFrec," | [",nXI-1,"] | "
      end subroutine subdivideTopo
!     subroutine allocintegPoints(iJ)
!     use soilVars, only : BETA
!     use resultVars, only : allpoints,BouPoints,npts,nBpts
!     use waveNumVars, only : DFREC
!     use Gquadrature, only: WLmulti!,Gqu_n => Gquad_n
!     implicit none
!     integer, intent(in) :: iJ
!     integer :: J,i_X,iXi
!     real :: MinWaveLenght,distance
!     J = iJ
!     if (J .lt. 2) J = 2
!     do iXi = 1,nBpts
!     do i_x = 1,nBpts
!     !Si la distancia entre segmentos es menor a CERCA usamos la quadratura
!       MinWaveLenght = min(BETA(BouPoints(i_X)%layer),BETA(BouPoints(iXi)%layer)) / (DFREC*real(J-1))
!       distance = sqrt((BouPoints(i_X)%center(1) - BouPoints(iXi)%center(1))**2 + &
!                       (BouPoints(i_X)%center(2) - BouPoints(iXi)%center(2))**2)
!       if (distance .lt. MinWaveLenght * WLmulti) then !use quadratura
!        BouPoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .true.
!      else !use integración directa
!        BouPoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .false.
!       end if
!     end do !i_x
!     do i_x = 1,npts
!     !Si la distancia entre segmentos es menor a CERCA usamos la quadratura
!       MinWaveLenght = min(BETA(allpoints(i_X)%layer),BETA(BouPoints(iXi)%layer)) / (DFREC*real(J-1))
!       distance = sqrt((allpoints(i_X)%center(1) - BouPoints(iXi)%center(1))**2 + &
!                       (allpoints(i_X)%center(2) - BouPoints(iXi)%center(2))**2)
!       if (distance .lt. MinWaveLenght * WLmulti) then !use quadratura
!        allpoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .true.
!      else !use integración directa
!        allpoints(i_X)%GT_k(iXi)%shouldUseQuadrature = .false.
!       end if
!     end do !i_x
!     end do !ixi
!     end subroutine allocintegPoints
      function interpol(xA,yA,xB,yB,x)
      implicit none
      real :: interpol
      real, intent(in) :: xA,yA,xB,yB,x
      real :: m
      !interpolación lineal de la ordenada de x que está entre A y B
      m = (yB-yA) / (xB-xA)
      interpol = yA + m * (x-xA)
      end  function interpol
      function thereisavirtualsourceat(iz)
      use resultVars, only : pota
      integer :: iz
      logical :: thereisavirtualsourceat
      thereisavirtualsourceat = .false.
      if (pota(iz,2) .gt. 0) then
         thereisavirtualsourceat = .true.
      end if
      end function thereisavirtualsourceat
      subroutine waveAmplitude(outpf)                                   !    FIX NFREC
      use wavelets !las funciones: ricker, fork
      use waveVars, only : Dt,Uo, tipoPulso
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
      end if!
      if (tipoPulso .eq. 0) then
        allocate(Uo(NFREC*2))
        Uo=cmplx(1,0,8)
      else
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
       call plotXYcomp(Uo,real(Dt,4),size(Uo),titleN,xAx,yAx,1200,800)
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
       call plotSpectrum(Uo,real(DFREC,4),size(Uo),int(size(Uo)/2),titleN,xAx,yAx,logflag,1200,800)
       CALL chdir("..")
       print*,"plotted "
!      logflag = 'logy     '
!      call plotSpectrum(Uo,DFREC,size(Uo),int(size(Uo)/2),titleN,xAx,yAx,logflag,1200,800)
                        
!      call SYSTEM('../plotSpectrum rick_frec.txt rick_frec.pdf frecuency[Hz] amplitude logx 1200 800')
!      !CALL SYSTEM ('../plotXYcomp rick_frec.txt & 
!      !            rick_frec_ugly.pdf time[sec] amplitude 1200 800')
      end if
      end if !tipo
      
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
      subroutine crunch(i_zfuente,J,cOME,outpf)
      ! esta función es llamada con cada profundidad donde hay
      ! por lo menos una fuente.
      use gloVars, only: verbose,PI,plotStress
      use resultVars, only : pota,Punto,nZs, & 
                             MecaElem,Hf,Hk, NxXxMpts, & 
                             trac0vec,ibemMat, NxXxMpts,FFres
      use Gquadrature, only : Gquad_n
      use refSolMatrixVars, only : B,Ak
      use waveNumVars, only : NMAX,DK, minKvalueW,minKvalueU,lapse,KtrimStart
      use wavelets !fork
      use meshVars, only: MeshDXmultiplo
      use meshVars, only: DeltaX, MeshDXmultiplo
      use dislin
      implicit none
      
      interface !porque 'asociar' recibe un puntero como argumento
        subroutine asociar(px, itabla_z, itabla_x)
          use resultVars, only : Punto
          type(Punto), pointer :: PX
          integer :: itabla_x, itabla_z
        end subroutine asociar
        
        subroutine doesanyoneneedsGQ(anyone,zi,ei,nPXs, & 
                                   itabla_z,& 
                                   xf,zf,ef)
          use resultVars, only : Punto
          logical, intent(out) :: anyone
          real, intent(out) :: zi
          integer, intent(out) :: ei
!         type(Punto), pointer :: PXs(:)
          integer, intent(out) :: nPXs
          integer, intent(in) :: itabla_z
!         integer, intent(in) :: nx
          real, intent(in) :: xf,zf
          integer, intent(in) :: ef
        end subroutine doesanyoneneedsGQ
        
        function Traction(RW,normal,l)
          complex*16 :: Traction,RW(:)
          real*8, intent(in), dimension(2) :: normal
          integer, intent(in) :: l
        end function Traction
      end interface
      
      
      integer, intent(in) :: i_zfuente,J
      complex*16, intent(in)  :: cOME
      integer, intent(in) :: outpf
      integer :: iGq,iGq_med,mecStart,mecEnd,xXx,dirStart,dirEnd
      integer :: itabla_z,itabla_x,ik,ef,dir,nGqs,ei,iMec,nPXs,ixi,nxis,i,thislapse
      real :: k,zf,xf,zi,factor,peso,mov_x
      logical :: intf
      type(Punto), pointer :: pXi,p_x
!     type(Punto), pointer :: pXs(:)
      logical :: anyone,warning
      real*8, dimension(2) :: nf
      type(MecaElem)  :: calcMecaAt_k_zi, Meca_diff!, Meca_ff,calcff
      complex*16, dimension(2*nmax,5), target :: auxK,savedAuxK
      integer, dimension(5,2) :: sign
      complex*16, dimension (:), pointer :: RW
      type(FFres),target :: FF
!     logical :: exist
!     character(LEN=90) :: nombre
      !  fza horiz     fza vert       !(para la crepa)
      sign(1,1) = -1 ; sign(1,2) = 1  !W        impar           par
      sign(2,1) = 1  ; sign(2,2) = -1 !U          par           impar
      sign(3,1) = -1 ; sign(3,2) = 1  !s33      impar           par
      sign(4,1) = 1  ; sign(4,2) = -1 !s31        par           impar
      sign(5,1) = -1 ; sign(5,2) = 1  !s11      impar           par
      factor = sqrt(2.*NMAX*2)
      
      ! indice itabla_x de la fuente en la tabla
      if (i_zfuente .eq. 0) then; itabla_x = 3 ! la fuente real
      else; itabla_x = 2 + pota(i_zfuente,1) + 1 ! la primer fuente virtual
      end if
      call asociar(px = pXi, itabla_z  = i_zfuente, itabla_x  = itabla_x)
      ef = pXi%layer
      intf = pXi%isOnInterface
      nf(1) = pXi%normal(1)
      nf(2) = pXi%normal(2)
      
      if (i_zfuente .eq. 0) then; nGqs = 1
      else; nGqs = Gquad_n
      end if
      iGq_med = ceiling(nGqs/2.0) ! nGqs +1   i de Punto para integración lineal.
      
      dirStart = 1; dirEnd = 2
      if (i_zfuente .eq. 0) then !ahorramos si la fuente tiene un componente nulo
        if (abs(nf(1)) .lt. 0.001) then
          dirStart = 2
        end if!
        if (abs(nf(2)) .lt. 0.001) then
          dirEnd = 1
        end if
      end if!      
      if (verbose .ge. 3) print*,"dirs", dirStart," -> ",dirEnd
      

      do dir =dirStart,dirEnd! componete de Fuerza horizontal; vertical..........
      do iGq = 1,nGqs ! cada punto de integración / centro ......................
        if (verbose .ge. 3) print*,"iGq ",iGq,"/",nGqs
      ! amplitud de las ondas planas dada la profundidad de la fuente
        if (i_zfuente .eq. 0) then
         xf = pXi%center(1) 
         zf = pXi%center(2)
        else
         xf = pXi%Gq_xXx_coords(iGq,1) ! = xfsource si iz = 0
         zf = pXi%Gq_xXx_coords(iGq,2) ! = zfsource si iz = 0
        end if!
        if (verbose .ge. 3) print*,"xi : (",xf,",",zf,")"
        
       do ik = 1,nmax+1 ! amplitud de ondas planas dada fuente en Z .......
        k = real(ik-1) * dK; if (ik .eq. 1) k = dk * 0.01
        call vectorB_force(B(:,ik),zf,ef,intf,dir,cOME,k)  
        B(:,ik) = matmul(AK(:,:,ik),B(:,ik))
       end do ! ik
      
      ! resultados para cada Z donde hay receptores ..............................
      do itabla_z = 1,nZs !(en la tabla pota todos los puntos son receptores)        
      ! ¿En este Z hay algún punto que requiera integración gaussiana?
        ! si : acumulamos resultado con este iGq
        ! no : esperamos a que iGq sea iGq_med
      call doesanyoneneedsGQ(anyone,zi,ei,nPXs, & 
                               itabla_z, & 
                               xf,zf,ef)
        if (verbose .ge. 3) write(outpf,*) "receptores en itabla_z = ", & 
                                           itabla_z, " -> ",zi,"m"
            ! con este if logramos excluir calculos innecesarios cuando
            ! a esta profundidad sólo existen receptores que requieren
            ! integración lineal dada una fuente virtual.
        if (anyone .or. (iGq .eq. iGq_med)) then
            ! ¿calcular sólo G, sólo T o ambos?
            mecStart = 1; mecEnd = 5
            if(pota(itabla_z,2) .eq. 0) mecEnd = 2 ! no T
            if(pota(itabla_z,1) .eq. 0) mecStart = 3 ! no G
            if (plotStress .and. (i_zfuente .eq. 0)) then
              mecStart = 1; mecEnd = 5
            end if
            if (verbose .ge. 3) print*,"mec (",mecstart,mecend,")"
            auxK = cmplx(0.0d0,0.0d0,8); savedauxk = auxK
            warning = .false.
            thislapse = lapse
            
          do ik = 1,nmax+1 ! k loop : campo difractado por estratos
            k = real(ik-1) * dK; if (ik .eq. 1) k = dk * 0.01
            Meca_diff = calcMecaAt_k_zi(B(:,ik),&  
                     zi,ei,cOME,k, & 
                     dir,mecStart,mecEnd,outpf)
            auxK(ik,mecStart:mecEnd) = Meca_diff%Rw(mecStart:mecEnd) 
            
            cycle
            ! trim K integration domain if it is zero 
            if (abs(auxK(ik,mecStart)) .lt. minKvalueW .and. & 
                abs(auxK(ik,mecEnd)) .lt. minKvalueU .and. & 
                warning .eqv. .true.) then 
                   if (verbose .ge. 1) print*,"trim loop at ik = ",ik
                   KtrimStart = max(KtrimStart,ik)
                   exit ! ----------------------------- esto es arriesgado 
            end if!
            if (abs(auxK(ik,mecStart)) .lt. minKvalueW .and. & 
                abs(auxK(ik,mecEnd)) .lt. minKvalueU) then 
                   thislapse = thislapse - 1
            else
                   if (thislapse .lt. lapse) &
                   thislapse = thislapse + 1
            end if!
            if (thislapse .le. 0) warning = .true.
          end do ! ik
          
         
        
      ! doblar la crepa ...............................................
           !        calculado              crepa 
           !   _________|____________ _______|________  
           !    1 2 3 ... NMAX NMAX+1  ...      2*NMAX     
           ! k=                       nmax*  ... (2)*
        do iMec = mecStart,mecEnd
         auxK(nmax+2:nmax*2,iMec) = sign(iMec,dir)* auxK(nmax:2:-1,iMec)         
        end do !iMec
        savedAuxK(1:2*nmax,mecStart:mecEnd) = auxK(1:2*nmax,mecStart:mecEnd)
        
        ! Coordenada X ....................
        if (i_zfuente .eq. 0) then
          nXis = 1 ! la fuente real es una
        else
          nXis = pota(i_zfuente,2) ! puede haber varias fuentes virtuales
        end if
        ! para cada fuente a la profundidad iz ..................................
       
        do iXi = 1,nXis
         if (i_zfuente .ne. 0) then ! si no es la fuente real, son las virtuales
             itabla_x = 2 + pota(i_zfuente,1) + iXi
             call asociar(px = pXi, itabla_z  = i_zfuente, itabla_x  = itabla_x)
             xf = pXi%Gq_xXx_coords(iGq,1)
             nf(dir) = pXi%normal(dir)
         end if!
         
        ! para cada X receptor a la profundidad itabla_z ......................
         if (verbose .ge. 3) print*,"nPXs ",nPXs
            do itabla_x =1,nPXs
              ! reponer auxK (cambia cada ciclo)
              auxK = savedAuxK
              call asociar(px = p_x, itabla_z  = itabla_z, itabla_x  = 2+ itabla_x)
              
              if (p_x%isBoundary) then
                if (pXi%boundaryIndex .eq. p_x%boundaryIndex)  then
                ! fuente y receptor son virtules y están en el mismo lugar
                    if (dir .eq. 1) then
                       ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                                   pXi%boundaryIndex *2 -(2 - dir)) = 0.5
                       ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                                   pXi%boundaryIndex *2 -(2 - dir)) = 0            
                    else ! dir = 2
                       ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                                   pXi%boundaryIndex *2 -(2 - dir)) = 0
                       ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                                   pXi%boundaryIndex *2 -(2 - dir)) = 0.5
                    end if
                    cycle !skip to the next receptor if any
                end if
              end if
                          
              ! información X de fuente y receptor
              do imec = mecStart,mecEnd
                auxk(1:nmax,imec) = auxk(1:nmax,imec) * & 
                exp(cmplx(0.0,(/((ik-1)*dk,ik=1,nmax)/)* & 
                     (p_x%center(1)-xf) ,8))
                     
                auxk(1     ,imec) = auxk(1,imec)* &
                exp(cmplx(0.0,0.01*dk * & 
                    (p_x%center(1)-xf) ,8))
        
                auxk(nmax+1:nmax*2,imec) = auxk(nmax+1:nmax*2,imec) * & 
                exp(cmplx(0.0,(/((-ik)*dk,ik=nmax,1,-1)/)* & 
                    (p_x%center(1)-xf) ,8))

                auxk(:,iMec) = auxk(:,iMec)* & 
                    cmplx(Hf(J)* Hk,0.0,8)/(4*pi)! taper 
              end do !imec

                      
!          inquire(file="out1.txt", exist=exist)
!          if (exist) then
!      open(12, file="out1.txt", status="old", position="append", action="write")
!          else
!      open(12, file="out1.txt", status="new", action="write")
!          end if
!          do ik=1,size(auxK(:,1))
!          write(12,'(E12.3,1x)',advance='NO') abs(auxK(ik,2))
!          end do !ik
!          write(12,*) " "
!          close(12)  
                     
              ! guardar el FK si es que es interesante verlo
              if ((mecStart .eq. 1) .and. (mecEnd .ge. 2)) then
              if(i_zfuente .eq. 0) then ; if (p_x%guardarFK) then
                if (dir .eq. dirStart) then
                  p_x%FK(J,1:nmax,1) = 0
                  p_x%FK(J,1:nmax,2) = 0
                end if
                
                p_x%FK(J,1:nmax,1) = &
                p_x%FK(J,1:nmax,1) + auxK(1:nmax,1)* nf(dir)
                p_x%FK(J,1:nmax,2) = &
                p_x%FK(J,1:nmax,2) + auxK(1:nmax,2)* nf(dir)
              end if ; end if !guardar
              end if
              
              ! K -> X
              auxK = auxK * factor
              do iMec = mecStart,mecEnd
                call fork(2*nmax,auxK(:,iMec),+1,verbose,outpf)
              end do
              auxK = auxK / factor
                 
              !          ,--- k  ->  x = 0 (x del receptor)
              !          |
              RW => auxK(1,1:5)
              
      ! Repartir resultado         
              ! decidimos que hacer con el resultado dependiendo del
              ! tipo de receptor y el tipo de fuente
              !
              !                    | allpoint  |   boupoint
              !   ----------------------------------------------
              !   fuente real      |     W,U   |  term. indep.
              !   fuente virtual   |      G    |      T
              
      if (i_zfuente .eq. 0) then !fuente real (no hay integral GQ)................
          if (p_x%isBoundary) then ! [term. indep.] ..............................
                     !  | Tx |
                     !  | Tz |
                     if (verbose .ge. 3) print*,"term indep"
                     ! el indice en el vector de terms indep depende
                     ! del punto receptor sobre la frontera
                     FF%W = 0;FF%U = 0;FF%Tx = 0;FF%Tz = 0
                     if (.not. pXi%isOnInterface) then
                     if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                        call calcFreeField(FF,dir,p_X%center,p_X%normal, & 
                        pXi%center,p_x%layer,cOME,mecStart,mecEnd)
                     end if
                     end if
                     trac0vec(p_x%boundaryIndex * 2 - (1 - 0)) = &
                   trac0vec(p_x%boundaryIndex * 2 - (1 - 0)) - &
                    (Traction(RW * nf(dir), p_x%normal,0) + FF%Tx * nf(dir))
                   
                   trac0vec(p_x%boundaryIndex * 2 - (1 - 1)) = &
                   trac0vec(p_x%boundaryIndex * 2 - (1 - 1)) - &
                    (Traction(RW * nf(dir), p_x%normal,1) + FF%Tz * nf(dir))
          else ! not boundary point [W,U] ........................................
            if (p_x%guardarMovieSiblings) then
!                      if (verbose .ge. 3) print*,"movie"
                       ! (1) reordenar la crepa existente en X
                     do iMec = mecStart,mecEnd                                    ! 1,5 -- 2
                      auxK(1:nmax*2,iMec) = cshift(auxK(1:nmax*2,iMec),SHIFT=nmax)
                     end do
                       ! (2) guardar sólo los interesantes
                       xXx = 1
                       do i=nmax+1-(nxXxmpts-1)/2* meshdxmultiplo, & 
                            nmax+1+(nxXxmpts-1)/2* meshdxmultiplo, & 
                            meshdxmultiplo
                            FF%W = 0;FF%U = 0;FF%Tz = 0;FF%Tx = 0
                            if (.not. pXi%isOnInterface) then
                            if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                              ! la coord X de este hermanito
                              mov_x = (i-(nmax+1)) * DeltaX 
                              if (verbose .ge. 3) print*,xXx,i,"--",mov_x
                              call calcFreeField(FF,dir,(/mov_x,p_x%center(2)/), & 
                             p_X%normal,pXi%center,p_x%layer,cOME,mecStart,mecEnd)
                            end if
                            end if
                            !FF%W = 0;FF%U=0
              p_x%WmovieSiblings(xXx,J,1) = (auxK(i,1) + FF%W)* nf(dir) + &
              p_x%WmovieSiblings(xXx,J,1)
              
              p_x%WmovieSiblings(xXx,J,2) = (auxK(i,2) + FF%U)* nf(dir) + &
              p_x%WmovieSiblings(xXx,J,2)              
            if (plotStress) then
              ! tZ
              p_x%WmovieSiblings(xXx,J,3) =  & 
              p_x%WmovieSiblings(xXx,J,3) + (auxK(i,4)*p_X%normal(1) + & 
                                         auxK(i,3)*p_X%normal(2)+ FF%Tz)* nf(dir)
              ! tX
              p_x%WmovieSiblings(xXx,J,4) =  & 
              p_x%WmovieSiblings(xXx,J,4) + (auxK(i,5)*p_X%normal(1) + & 
                                         auxK(i,4)*p_X%normal(2)+ FF%Tx)* nf(dir)
            end if        
                            xXx = xXx + 1
                       end do !i de pelicula
           else !not a movie point ...............................................
                       FF%W = 0;FF%U = 0;FF%Tz = 0;FF%Tx = 0
                     if (.not. pXi%isOnInterface) then
                     if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                     call calcFreeField(FF,dir,p_X%center,p_X%normal, & 
                        pXi%center,p_x%layer,cOME,mecStart,mecEnd) 
                     end if
                     end if                                                 !
                       p_x%W(J,1) = &
                       p_x%W(J,1) + (RW(1)+ FF%W) * nf(dir)
                       p_x%W(J,2) = &
                       p_x%W(J,2) + (RW(2)+ FF%U) * nf(dir)
                       
               if (plotStress) then
               p_x%W(J,3) =  & 
               p_x%W(J,3) + (RW(4)*p_X%normal(1) + & 
               RW(3)*p_X%normal(2)+ FF%Tz)* nf(dir)
               p_x%W(J,4) =  & 
               p_x%W(J,4) + (RW(5)*p_X%normal(1) + & 
               RW(4)*p_X%normal(2)+ FF%Tx)* nf(dir)
               end if
           end if ! movie point or not
        end if !boundary point or not
      else !fuente virtual (puede haber integral GQ)..............................
                if (nGqs .eq. 1) then
                  peso = pXi%length
                else
                  peso = pXi%Gq_xXx_C(iGq)
                end if!
        if (p_x%isBoundary) then ! [T] ...........................................
                     ! obtener func. Green de tracciones y
                     ! llenar la matriz del ibem
                     
                     !  |  Txx Txz  |   
                     !  |  Tzx Tzz  |
                     ! macro columna: pX%boundaryIndex         (fuente)
                     ! miniCol::dir 1 -> X ; 2 -> Z          (dir fuente)
                     ! macro renglón: p_x%boundaryIndex   (receptor)
                     ! miniRow::                 (componente de tracción)
                     FF%Tx = 0; FF%Tz = 0
                    if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                       call calcFreeField(FF,dir,p_X%center,p_X%normal, & 
                        (/xf,zf/),p_x%layer,cOME,mecStart,mecEnd)
                     end if
                     
                ! tracción en el receptor con la normal del receptor
                ! multiplicado por el peso de integración de la fuente
                     ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                             pXi%boundaryIndex *2 -(2 - dir)) = &
                     ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                             pXi%boundaryIndex *2 -(2 - dir)) + &
                       ((Traction(RW, p_x%normal,0) + FF%Tx) * peso)
                       
                    ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                            pXi%boundaryIndex *2 -(2 - dir)) = &
                    ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                            pXi%boundaryIndex *2 -(2 - dir)) + &
                       ((Traction(RW, p_x%normal,1) + FF%Tz) * peso)
                     
        else ! [G] ...............................................................
            if (p_x%guardarMovieSiblings) then ! .................................
                       ! (1) reordenar la crepa existente en X
                       do iMec = 1,2
                        auxK(1:nmax*2,iMec) = cshift(auxK(1:nmax*2,iMec),SHIFT=nmax)
                       end do
                       ! (2) guardar sólo los interesantes
                       xXx = 1
                       do i=nmax+1-(nxXxmpts-1)/2* meshdxmultiplo, & 
                            nmax+1+(nxXxmpts-1)/2* meshdxmultiplo, & 
                            meshdxmultiplo
                            FF%W = 0
                            FF%U = 0
                            if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                              ! la coord X de este hermanito
                              mov_x = (i-(nmax+1)) * DeltaX 
                              if (verbose .ge. 3) print*,i,"--",mov_x
                              call calcFreeField(FF,dir,(/mov_x,p_x%center(2)/),p_X%normal, & 
                              (/xf,zf/),p_x%layer,cOME,mecStart,mecEnd)
                            end if
                       ! func de green (sólo desplazamientos)
                            p_x%GT_gq_mov(pXi%boundaryIndex,1,dir,xXx) = & 
                            p_x%GT_gq_mov(pXi%boundaryIndex,1,dir,xXx) + &
                            (auxK(i,1) + FF%W) * peso
                            p_x%GT_gq_mov(pXi%boundaryIndex,2,dir,xXx) = & 
                            p_x%GT_gq_mov(pXi%boundaryIndex,2,dir,xXx) + &
                            (auxK(i,2) + FF%U) * peso
            
                            xXx = xXx + 1
                       end do
             else !not a movie point .............................................
                        FF%W = 0
                        FF%U = 0
                     if (p_x%layer .eq. pXi%layer) then !agregar campo libre
                     call calcFreeField(FF,dir,p_x%center,p_X%normal, & 
                              (/xf,zf/),p_x%layer,cOME,mecStart,mecEnd)
                     end if
                       p_x%GT_gq(pXi%boundaryIndex,1,dir) = &
                       p_x%GT_gq(pXi%boundaryIndex,1,dir) + &
                       (RW(1) + FF%W) * peso
                       p_x%GT_gq(pXi%boundaryIndex,2,dir) = &
                       p_x%GT_gq(pXi%boundaryIndex,2,dir) + &
                       (RW(2) + FF%U) * peso
                ! desplazamientos en el receptor
                ! multiplicado por el peso de integración de la fuente
             end if ! movie point or not?
         end if ! is it boundary or not?
      end if !tipo de fuente (real o virtual)
      end do ! itabla_x (receptor/columna)
      end do !iXi (fuente/renglón)
          
      end if ! anyone      
      end do ! itabla_z (receptor)
      
      
      end do ! iGq (fuente)
!     stop "crunch"
      end do ! dir
      
      
      
      end subroutine crunch
      
      subroutine doesAnyOneNeedsGQ(anyone,zi,ei,nPXs, & 
                                   itabla_z, & 
                                   xf,zf,ef)
                                   
      use resultvars, only : pota,punto,nBpts,npts
      use soilVars, only : beta
      use waveNumVars, only : Frec
      use Gquadrature, only: WLmulti
!     use glovars, only : workboundary
      implicit none
      
      interface !porque 'asociar' recibe un puntero como argumento
        subroutine asociar(px, itabla_z, itabla_x)
          use resultVars, only : Punto
          type(Punto), pointer :: PX
          integer :: itabla_x, itabla_z
        end subroutine asociar
      end interface
      
      type(Punto), pointer :: PX
!     type(Punto), pointer :: PXs(:)
      logical, intent(out) :: anyone
      integer, intent(out) :: ei,nPXs
      integer, intent(in) :: itabla_z,ef
      real, intent(out) :: zi
      integer :: ix,maxix
      real, intent(in) :: xf,zf
      real*8 :: MinWaveLenght,distance
      
!     needsGQ = .false.
      anyone = .false.
!     if (associated(PXs)) then 
!        nullify(PXs)
!     end if!
      if (itabla_z .eq. 0) then 
         call asociar(px,itabla_z,0)
         zi = PX%center(2)
         ei = PX%layer
         nPXs = nPts + nBpts
!        allocate(pXs(npxs))
!        do ix = 1,nPts
!          pxs(1 : npts) => allpoints(:)
!        end do
!        if (workboundary) then
!          do ix = 1,nBpts
!            pxs(1 + npts : nbpts + npts) => boupoints(:)
!          end do
!        end if
         
!        if (workboundary) then
!           PXs = (/ allpoints(:),boupoints(:) /)
!        else
!           PXs => allpoints(:)
!        end if
        
         return
      end if 
      maxix = pota(itabla_z,1) + pota(itabla_z,2)
!     allocate(PXs(maxix)) !sin las dos columnas iniciales de indices
      nPXs = maxix
      do iX = 3,maxix+2
        call asociar(PX,itabla_z,ix)
!       PXs(ix-2) = PX
        MinWaveLenght = min(real(abs(BETA(PX%layer)),4),real(abs(BETA(ef)),4)) / FREC
        distance = sqrt((PX%center(1) - xf)**2 + &
                        (PX%center(2) - zf)**2)
        if (distance .lt. MinWaveLenght * WLmulti) then !use quadratura
          anyone = .true.
          exit
        end if
      end do !ix
      zi = PX%center(2)
      ei = PX%layer
!     print*,"anyone:", needsGQ
      end subroutine doesanyoneneedsGQ
      
      subroutine asociar(PX,itabla_z, itabla_x)
      use resultVars, only : pota,Punto,allpoints,boupoints,Po
      type(Punto), pointer :: PX
      integer :: itabla_x, itabla_z
      ! apuntador al punto en la tabla
!       if (associated(PX)) then 
          nullify(PX)
!       end if!
        if (itabla_z .ne. 0) then
          if (pota(itabla_z, itabla_x) .gt. 0) then
             PX => allpoints(pota(itabla_z, itabla_x))!; print*,"allp "
          else
             PX => boupoints(abs(pota(itabla_z, itabla_x)))!; print*,"boup "
          end if
        else
          PX => Po
        end if
!      print*,"sucess asociar"
      end subroutine asociar
      subroutine matrixA_borderCond(this_A,k,cOME_i,outpf)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : verbose,UI,UR,Z0!,PI!,dp
!     use waveVars, only : Theta     
!     use refSolMatrixVars, only : A ! aquí guardamos el resultado 
      use debugStuff  
!     use waveNumVars, only : K !DWN compliant / explicitly
      implicit none
      
!     real, parameter    :: x_i = 0
      complex*16,    intent(inout), dimension(4*N+2,4*N+2) :: this_A
      real*8               :: z_i
      real*8,       intent(in)     :: k
      complex*16, intent(in)     :: cOME_i  
      
      integer, intent(in) :: outpf

      complex*16, dimension(2,4) :: subMatD, subMatS, subMatD0, subMatS0
!     complex*16, dimension(4,1) :: diagMat
      complex*16, dimension(4,4) :: diagMat
      complex*16 :: gamma,nu,xi!,mukanu2,mukagam2,l2m,g2,k2,lk2!,lg2!,eta
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
          
!         mukanu2  = 2* amu(e)* k * nu
!         mukagam2 = 2* amu(e)* k * gamma
!         l2m = lambda(e) + 2 * amu(e)
          xi = k**2-nu**2 !(nu**2 - k**2) * amu(e)
!         g2 = gamma**2
!         k2 = k**2
!         lk2 = lambda(e)*k2
!         lg2 = lambda(e)*g2
!         eta = 2*gamma**2 - cOME_i**2/BETA(e)**2
          

          ! en fortran los elementos se indican por columnas:
          subMatD0 = RESHAPE((/ -gamma,-k*UR,-k*UR,nu,& 
                                 gamma,-k*UR,-k*UR,-nu /), &
                           (/ 2,4 /))
          subMatD0 = UI * subMatD0 
          subMatS0 = RESHAPE((/ xi,-2*k*gamma,-2*k*nu,-xi,& 
                               xi,2*k*gamma,2*k*nu,-xi /),&
                           (/2,4/)) 
          subMatS0 = amu(e) * subMatS0                  
!         subMatS0 = RESHAPE((/ l2m*(-1*g2)-lk2 , -mukagam2, &
!                              -mukanu2        , xi, &
!                              l2m*(-1*g2)-lk2 , mukagam2, & 
!                              mukanu2         , xi /),&
!                          (/2,4/))
          
          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1 
        do bord = 0,1
          if (e+bord > N+1) then ! si 1+0;1+1;2+0;[2+1] > 2
            exit
          end if                           
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
!           diagMat = RESHAPE((/ egammaN, enuN, egammaP, enuP /), &
!                          (/ 4,1 /))
         diagMat = RESHAPE((/ egammaN, Z0, Z0, Z0, & 
                              Z0,    enuN, Z0, Z0, & 
                              Z0, Z0, egammaP, Z0, & 
                              Z0, Z0, Z0, enuP /), &
                           (/ 4,4 /))
            
          ! desplazamientos estrato i (en Fortran se llena por columnas)
!         subMatD = UI * subMatD !* exp(-UI * k * x_i) 
!         subMatD(:,1) = subMatD(:,1) * diagMat(1,1) !Down P
!         subMatD(:,2) = subMatD(:,2) * diagMat(2,1) !Down S
!         subMatD(:,3) = subMatD(:,3) * diagMat(3,1) !Up   P
!         subMatD(:,4) = subMatD(:,4) * diagMat(4,1) !Up   S
          
          subMatD = matmul(subMatD0,diagMat)
          
          ! esfuerzos estrato i
!         subMatS = AMU(e) * subMatS !* exp(-UI*k*x_i) 
!         subMatS(:,1) = subMatS(:,1) * diagMat(1,1) !Down P
!         subMatS(:,2) = subMatS(:,2) * diagMat(2,1) !Down S
!         subMatS(:,3) = subMatS(:,3) * diagMat(3,1) !Up   P
!         subMatS(:,4) = subMatS(:,4) * diagMat(4,1) !Up   S
          
          subMatS = matmul(subMatS0,diagMat)
          
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
             exit
           end if
          end if
          
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
            this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatD
            this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatS
          end if
          
        end do !bord loop del borde i superior o nferior
          iR= iR+4 
          iC= iC+4
          
      END DO !{e} loop de las macro columnas para cada estrato
          
          if (verbose >= 3) then
             call showMNmatrixZ(4*N+2,4*N+2,this_A,"A    ",outpf) 
          end if
      
      end subroutine matrixA_borderCond
      
      
      subroutine vectorB_force(this_B,z_f,e,fisInterf,direction,cOME,k)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO
      use gloVars, only : UR,UI,PI  
      use debugStuff   
!     use resultVars, only : allpoints,NPts
      implicit none
      
      complex*16,    intent(inout), dimension(4*N+2) :: this_B
      integer,    intent(in)    :: direction,e
      real,       intent(in)    :: k,z_f
      complex*16, intent(in)    :: cOME
      logical,    intent(in)    :: fisInterf
      integer :: iIf,nInterf!,iP
      real    :: SGNz!,x_i
      complex*16 :: gamma,nu,DEN,L2M
      real     :: errT = 0.001
      complex*16 :: omeAlf,omeBet!,AUX
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
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
          ! nencial a medida que z crece.
          if(aimag(gamma).gt.0.0)gamma = -gamma
          if(aimag(nu).gt.0.0)nu=-nu
          
      ! en cada interfaz (1 arriba) y (2 abajo)
      do iIf = 1,2
!         print*,"[",x,",",z_loc(iIf),"]"
!         if (z_loc(iIf) .ne. 0.) then
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
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
        s331(1) = 0
        S131(1) = - cmplx(1.0 / (2.0 * PI ),0.0,8)
      elseif (direction .eq. 2) then
        s333(1) = - cmplx(1.0 / (2.0 * PI ),0.0,8)
        S313(1) = 0
      end if
        G33(1) = 0
        G31(1) = 0
        G11(1) = 0
      end if
      
      ! El vector de términos independientes genera el campo difractado
      
      if (direction .eq. 1) then ! fuerza HORIZONTAL
      !                     =      (1) interfaz de arriba
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2) = + G31(1)!  w
        this_B(1+4*(e-1)-1) = + G11(1)!  u
       end if 
        this_B(1+4*(e-1)  ) = + S331(1)! szz
        this_B(1+4*(e-1)+1) = + S131(1)! szx   ! delta
    
      if (.not. fisInterf) then ! la fuerza no calló en la interfaz
      !                     =      (2) interfaz de abajo
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2) = - G31(2)!  w
        this_B(1+4*(e-1)+3) = - G11(2)!  u
        this_B(1+4*(e-1)+4) = - S331(2)! szz
        this_B(1+4*(e-1)+5) = - S131(2)! szx
       end if
      end if
      elseif (direction .eq. 2) then ! fuerza VERTICAL
      !                     =     (1) interfaz de arriba
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2) = + G33(1)!  w 
        this_B(1+4*(e-1)-1) = + G31(1)!  u 
       end if 
        this_B(1+4*(e-1)  ) = + S333(1)! szz   ! delta
        this_B(1+4*(e-1)+1) = + S313(1)! szx 
    
      if (.not. fisInterf) then
      !                     =      (2) interfaz de abajo
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2) = - G33(2)!  w 
        this_B(1+4*(e-1)+3) = - G31(2)!  u
        this_B(1+4*(e-1)+4) = - S333(2)! szz 
        this_B(1+4*(e-1)+5) = - S313(2)! szx 
       end if
      end if
      end if ! direction
!     print*,""
!     print*,this_B; stop "B"
      end subroutine vectorB_force
      function calcMecaAt_k_zi(thisIP_B,z_i,e,cOME_i,k,dir,mecStart,mecEnd,outpf)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : verbose,UI,UR,Z0!,PI
!     use waveVars, only : Theta
      use resultVars, only : MecaElem
      implicit none
      type (MecaElem)              :: calcMecaAt_k_zi
      real, intent(in)             ::  z_i, k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e,outpf,dir,mecStart,mecEnd
      complex*16, dimension(4*N+2),intent(in) :: thisIP_B
      complex*16 :: gamma,nu,xi,eta!,mukanu2,mukagam2,l2m,g2,k2,lk2,lg2
      complex*16 :: egammaN,enuN,egammaP,enuP
      complex*16, dimension(2,4) :: subMatD
      complex*16, dimension(3,4) :: subMatS
      complex*16, dimension(4,4) :: diagMat
      complex*16, dimension(4,1) :: partOfXX
      complex*16, dimension(2,1) :: resD
      complex*16, dimension(3,1) :: resS
      
      if (verbose >= 3) then
       write(outpf,'(a,F7.3,a,F12.7,a,F10.2,a,F10.2,a,I0)') & 
                    "finding solution values at w:", & 
                    real(cOME_i),"k=",k," {",z_i,"} e=",e
      end if
      
       
      ! algunas valores constantes para todo el estrato
      gamma = sqrt(cOME_i**2. /ALFA(e)**2. - k**2.)
      nu = sqrt(cOME_i**2. /BETA(e)**2. - k**2.)
      
      if(aimag(gamma).gt.0.0)gamma = -gamma
      if(aimag(nu).gt.0.0)nu=-nu
      
!         mukanu2  = 2* amu(e)* k * nu
!         mukagam2 = 2* amu(e)* k * gamma
!         l2m = lambda(e) + 2 * amu(e)
          xi = k**2-nu**2!(nu**2 - k**2) * amu(e)
!         g2 = gamma**2
!         k2 = k**2
!         lk2 = lambda(e)*k2
!         lg2 = lambda(e)*g2
          eta = 2*gamma**2 - cOME_i**2/BETA(e)**2
          
      !downward waves
          egammaN = exp(-UI * gamma * (z_i-Z(e)))
          enuN = exp(-UI * nu * (z_i-Z(e)))
          !upward waves 
          if (e /= N+1) then !(radiation condition)
            egammaP = exp(UI * gamma *(z_i-Z(e+1)))
            enuP = exp(UI * nu * (z_i-Z(e+1)))
          else
            egammaP = 0.0d0
            enuP = 0.0d0
          end if
          !la matrix diagonal
          diagMat = RESHAPE((/ egammaN, Z0, Z0, Z0, & 
                              Z0,    enuN, Z0, Z0, & 
                              Z0, Z0, egammaP, Z0, & 
                              Z0, Z0, Z0, enuP /), &
                           (/ 4,4 /))
      !coeficientes de las ondas en el estrato
        if (e /= N+1) then
          partOfXX(1:4,1) = thisIP_B(4*(e-1)+1 : 4*(e-1)+4)
        else !( condición de radiación)
          partOfXX(1:2,1) = thisIP_B(4*(e-1)+1 : 4*(e-1)+2)
          partOfXX(3:4,1) = (/z0,z0/)
        end if!
        if (verbose>=3) then
          print*,"xx1 ",partOfXX(1,1)
          print*,"xx2 ",partOfXX(2,1)
          print*,"xx3 ",partOfXX(3,1)
          print*,"xx4 ",partOfXX(4,1)
        end if  
      
      ! desplazamientos
      if (mecStart .eq. 1)then
        subMatD = RESHAPE((/ -gamma,-k*UR,-k*UR,nu, & 
                            gamma,-k*UR,-k*UR,-nu /), &
                           (/ 2,4 /))
        subMatD = UI * subMatD
        subMatD = matmul(subMatD,diagMat)
        resD = matmul(subMatD,partOfXX)
        
!       calcMecaAt_k_zi%Rw(1) = resD(1,1) !W
!       calcMecaAt_k_zi%Rw(2) = resD(2,1) !U
        if (dir .eq. 2) then !vertical
          calcMecaAt_k_zi%Rw(1) = resD(1,1) !W
          calcMecaAt_k_zi%Rw(2) = - resD(2,1) !U 
        else ! dir eq. 1 !horizontal
          calcMecaAt_k_zi%Rw(1) = - resD(1,1) !W
          calcMecaAt_k_zi%Rw(2) = resD(2,1) !U 
        end if
      end if ! desplazamientos
      
      ! esfuerzos
      if (mecEnd .eq. 5) then
     
      subMatS = RESHAPE((/ xi,      -2.0*k*gamma,     eta,     &
                          -2.0*k*nu,     -xi,        2.0*k*nu,   &
                           xi,       2.0*k*gamma,     eta,     &
                           2.0*k*nu,     -xi,       -2.0*k*nu /),&
                           (/3,4/))      
     
      subMatS = amu(e) * subMatS
!       subMatS = RESHAPE((/ l2m*(-1*g2)-lk2 , -mukagam2, l2m*(-1*k2)-lg2, &
!                              -mukanu2        , xi   , mukanu2, &
!                              l2m*(-1*g2)-lk2 , mukagam2, l2m*(-1*k2)-lg2, & 
!                              mukanu2         , xi   , -mukanu2 /),&
!                          (/3,4/))
        subMatS = matmul(subMatS,diagMat)
        resS = matmul(subMatS,partOfXX)
        
!       calcMecaAt_k_zi%Rw(3) = resS(1,1) !s33
!       calcMecaAt_k_zi%Rw(4) = resS(2,1) !s31
!       calcMecaAt_k_zi%Rw(5) = resS(3,1) !s11 
        if (dir .eq. 2) then ! vertical
          calcMecaAt_k_zi%Rw(3) = resS(1,1) !s33
          calcMecaAt_k_zi%Rw(4) = - resS(2,1) !s31
          calcMecaAt_k_zi%Rw(5) = resS(3,1) !s11
        else ! dir eq. 1 !horizontal
          calcMecaAt_k_zi%Rw(3) = - resS(1,1) !s33
          calcMecaAt_k_zi%Rw(4) = resS(2,1) !s31
          calcMecaAt_k_zi%Rw(5) = - resS(3,1) !s11
        end if
        
        
      end if ! esfuerzos
          
      end function calcMecaAt_k_zi
      
      subroutine calcFreeField(FF,j,p_X,nX,pXi,e,cOME,mecStart,mecEnd)
      !                      |  |   |  |__ estrato de fuente y receptor
      !                      |  |   |__ coordenadsa de la fuente (2)
      !                               normal al receptor
      !                      |  |__ coordenadas del receptor (2)
      !                      |__ dirección de la fuerza (1;2) :: (X,Z)
      ! solucion de campo libre para fuente y receptor en el mismo estrato
      
      ! Kaussel, Fundamental solutions in elastodynamics... pag 38
      use soilVars ,only : alfa,beta,amu,lambda,rho
      use gloVars  
!     use waveNumVars, only : NFREC,NMAX!,DK
      use hank
      use resultvars, only : FFres
      implicit none
      real,       intent(in)     :: p_X(2),pXi(2)
      integer,    intent(in)     :: j,e,mecStart,mecEnd
      complex*16, intent(in)     :: cOME
!     complex*16, dimension(5)   :: freef
      real*8,     intent(in)     :: nX(2)
      type (FFres), intent(out) :: FF
      real*8 :: r,gamma(2)
      complex*16 :: A,B,C,Dqr,Dkr
      complex*16 :: omeP,omeS!,psi,Dpsi,chi,Dchi !preliminary
!     complex*16 :: dg1j1,dg1j3,dg3j1,dg3j3 !complenetary
      complex*16 :: H0s,H1s,H2s,H0p,H1p,H2p !Hankel 
!     complex*16 :: i_4, b_a_2 !auxiliar
      integer :: i,deltaij!,k
      
!     complex*16 :: omeP,omeS,psi,Dpsi,chi,Dchi !preliminary
!     freef = 0
      ! preliminary functions
      r = sqrt((p_x(1)-pXi(1))**2 + (p_x(2)-pXi(2))**2)
      gamma(1) = (p_X(1) - pXi(1)) / r ! gamma x
      gamma(2) = (p_X(2) - pXi(2)) / r ! gamma z
!     print*,r
!     print*,gamma
!     print*,alfa(e)
!     print*,beta(e);stop 0
      omeP = cOME * r / alfa(e)
      omeS = cOME * r / beta(e)
!     i_4 = cmplx(0.0,0.25,8) ! i/4
!     b_a_2 = (beta(e)/alfa(e))**2 ! (beta/alfa)**2
      
!     print*,"(",p_x(1),p_x(2),") - (",pXi(1),pXi(2),") - ",omeP," : ",omeS
      ! funcs de Hankel de segunda especie
      call hankels(omeP,H0p,H1p)
      H2p = -H0p + 2/omeP * H1p
      call hankels(omeS,H0s,H1s)
      H2s = -H0s + 2/omeS * H1s
      
!     print*,"j",j
!     nX = (/1.0/(2.0**0.5),1.0/(2.0**0.5)/)
      
      A = H0p/alfa(e)**2 + H0s/beta(e)**2
      B = H2p/alfa(e)**2 - H2s/beta(e)**2
      
      if (mecStart .eq. 1) then
      ! W
      i = 2
      FF%W = -UI/8.0/rho(e)*(A*deltaij(i,j)-(2*gamma(i)*gamma(j)-deltaij(i,j))*B)
!     print*,"W == ",FF%W, ":",abs(FF%W)
      
      ! U
      i = 1
      FF%U = -UI/8.0/rho(e)*(A*deltaij(i,j)-(2*gamma(i)*gamma(j)-deltaij(i,j))*B)
!     print*,"U == ",FF%U, ":",abs(FF%U)
      end if!
      
      if (mecEnd .eq. 5) then
      Dqr = omeP*H1p
      Dkr = omeS*H1s
      C = Dqr/alfa(e)**2 -Dkr/beta(e)**2
      ! TZ
      i = 2
      FF%Tz = amu(e)*UI /(2*rho(e)*r)*((B+(lambda(e)*Dqr)/(2*amu(e)*alfa(e)**2))*gamma(j)*nX(i) &
      + (B + Dkr/(2*beta(e)**2))*(gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))
!     print*," Tz == ",FF%Tz, ":",abs(FF%Tz)
      
      ! TX
      i = 1
      FF%Tx = amu(e)*UI/(2*rho(e)*r)*((B+(lambda(e)*Dqr)/(2*amu(e)*alfa(e)**2))*gamma(j)*nX(i) &
      + (B + Dkr/(2*beta(e)**2))*(gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))
!     print*," Tx == ",FF%Tx, ":",abs(FF%Tx)
      end if
      
!      Kausell
!     ! psi, Dpsi, chi, Dchi
!     psi = i_4 * (((H1s/omeS) - b_a_2 * (H1p/omeP)) - H0s)
!     chi = i_4 * (b_a_2 * H2p - H2s)
!     
!     if (mecStart .eq. 1) then
!     ! W =>  i = 2     g3j
!     freef(1) = 1/amu(e)* (psi * deltaij(2,j) + chi * gamma(2) * gamma(j))
!     ! U =>  i = 1     g1j
!     freef(2) = 1/amu(e)* (psi * deltaij(1,j) + chi * gamma(1) * gamma(j))
!     end if!
!     if (mecEnd .eq. 5) then
!     
!     
!     
!     ! complementary
!     Dpsi = chi/r + i_4/r * omeS * H1s
!     Dchi = i_4/r * (b_a_2 * omeP * H1p - omeS * H1s) - 2 * chi / r
!     ! d(gij)/dxk
!     i = 1 !componente
!     k = 1 !derivada
!     dg1j1 = 1/amu(e) * (gamma(k)*(Dpsi * deltaij(i,j) + & 
!                              (Dpsi - 2 * chi / r)*gamma(i)*gamma(j)) + & 
!                    chi/r * (deltaij(i,k)*gamma(j) + deltaij(j,k)*gamma(i)))
!     print*,"dg1j1",dg1j1
!     i = 1 !componente
!     k = 2 !derivada
!     dg1j3 = 1/amu(e) * (gamma(k)*(Dpsi * deltaij(i,j) + & 
!                              (Dpsi - 2 * chi / r)*gamma(i)*gamma(j)) + & 
!                    chi/r * (deltaij(i,k)*gamma(j) + deltaij(j,k)*gamma(i)))
!     print*,"dg1j3", dg1j3
!     i = 2 !componente
!     k = 1 !derivada
!     dg3j1 = 1/amu(e) * (gamma(k)*(Dpsi * deltaij(i,j) + & 
!                              (Dpsi - 2 * chi / r)*gamma(i)*gamma(j)) + & 
!                    chi/r * (deltaij(i,k)*gamma(j) + deltaij(j,k)*gamma(i)))   
!     print*,"dg3j1", dg3j1
!     i = 2 !componente
!     k = 2 !derivada
!     dg3j3 = 1/amu(e) * (gamma(k)*(Dpsi * deltaij(i,j) + & 
!                              (Dpsi - 2 * chi / r)*gamma(i)*gamma(j)) + & 
!                    chi/r * (deltaij(i,k)*gamma(j) + deltaij(j,k)*gamma(i)))            
!     ! szz
!     print*,"dg3j3", dg3j3
!     print*,"lambda ",lambda(e)
!     print*,"mu ",amu(e)
!     freef(3) = lambda(e) * (dg3j3 + dg1j1) + 2* amu(e) * dg3j3
!     ! szx
!     freef(4) = amu(e) * (dg3j1 + dg1j3)
!     ! sxx
!     freef(5) = lambda(e) * (dg3j3 + dg1j1) + 2* amu(e) * dg1j1
!     end if
      end subroutine calcFreeField
      
      function deltaij(i,j)
      integer :: deltaij,i,j
      ! i : dirección del receptor {1;2} x,z
      ! j : dirección de la fuente {1;2} x,z
      
      deltaij = 0
      ! 1,1
      if (i .eq. j) then
        deltaij = 1
!       return
      end if
!     print*,"delta",i,j,"=",deltaij
      ! 3,2
!     if ((i .eq. 3) .and. (j .eq. 2)) then
!       deltaij = 1
!       return
!     end if 
      end function deltaij
      
      function calcFF(z_f,e_f,dir,zX,eX,k,cOME)
      use soilVars 
      use gloVars  
!     use waveNumVars, only : NFREC,NMAX!,DK
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
      complex*16                 :: gamma,nu,DEN,L2M
      complex*16                 :: omeAlf,omeBet
      real     :: errT = 0.001
      egamz=cmplx(1.0,0.0,8);enuz=cmplx(1.0,0.0,8);SGNz=0
      calcFF%Rw = cmplx(0.0,0.0,8)
      
      if (e_f .ne. eX) return
      
      DEN = 4.0*PI*RHO(e_f)*cOME**2.0
      omeAlf = cOME**2/ALFA(e_f)**2.0
      omeBet = cOME**2/BETA(e_f)**2.0
      L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
      
      z_i = zX - z_f
!     if (z_i .ne. 0) then
      if (abs(z_i) > errT) then
        SGNz = z_i / ABS(z_i)
      end if
      
        gamma = sqrt(omeAlf - k**2.0)
        nu = sqrt(omeBet - k**2.0)
        if(aimag(gamma).gt.0.0) then 
           gamma = -gamma
        end if!
        if(aimag(nu).gt.0.0) then 
           nu=-nu
        end if 
        
        
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
     
        calcFF%RW(1) = G31 !W
        calcFF%RW(2) = G11 !U
        calcFF%RW(3) = S331!s33
        calcFF%RW(4) = S131!s31
        calcFF%RW(5) = S111!s11
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
      
!     call dpsv
      
      
      call zgetrf(n,n,A,n,ipiv,info)
      if(info .ne. 0) stop "Problem at LU factorization of matrix "
      
      
      call zgetri(n,a,n,ipiv,work,lwork,info)
      if(info .ne. 0) stop "Problem at inverse of matrix "
      
      deallocate(work)
      deallocate(ipiv)
      end subroutine inverseA
      subroutine makeTaperFuncs_cutF_cutK(Npol,cutoff_fracFmax,cutoff_fracKmax)
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
         
         Hk(1:nmax+1) = (/ (i,i=0,nmax) /) * DK
         Hk(nmax+2:2*nmax) = (/ (i,i=nmax-1,1,-1) /) * DK 
         HK = 1.0/sqrt(1.0+(Hk/(cutoff_fracKmax*NMAX*DK))**(2*Npol))
         
      end subroutine makeTaperFuncs_cutF_cutK  
      
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
      expK = exp(cmplx(0.0 ,-2*pi * aux,8))
      
!     call showMNmatrixZ (2*nmax,1,expK,"nmax ",6)
      end subroutine expPIK
      function Traction(RW,normal,l)
      implicit none
      complex*16 :: Traction,RW(:)
      real*8, intent(in), dimension(2) :: normal
      integer, intent(in) :: l
!     print*,"got this RW ",RW
!     print*,"normal = ",normal
      Traction = RW(5-l) * normal(1) + &
                 RW(4-l) * normal(2)
      ! las tracciones:
!     print*,"Traction: ",Traction
      !   ,--- componente de la tracción : m
      !   |     ,--- cara
      !   |     |
      !  Tx = Sxx nx + Szx nz  |___ (0) fza real (fuente)  .
      !  Tz = Szx nx + Szz nz  |                           .  
      !
      !   ,--- componente de la tracción : l
      !   |,--- (1),(3) dirección de la fuerza : m
      !   ||    ,--- cara
      !   ||    |,--- fza
      !  Txx = Sxx1 nx + Szx1 nz  |___ (1) fza horizontal  .
      !  Tzx = Szx1 nx + Szz1 nz  |                        . 
      !  Txz = Sxx3 nx + Szx3 nz |___ (3) fza vertical     .
      !  Tzz = Szx3 nx + Szz3 nz |                         .
      
      !  T_lm = S_lkm * n_k
      !       = S_l1m * n_1 + S_l3m * n3
      !         __|__         __|__
      !      s11     s31   s13    s33
      !      s11     s31   s31    s33  (son equivalentes)
      !       5       4     4      3   (indice en RW(_,i) )
      !       0       1     0      1   (indice de submatriz: l )
      end function Traction



         
 

!     function integralEq16(iP_x,l,iPxi,m)
!        !                              receptor---,       ,--- fuente
!        !                                        _|___   _|___
!        !  ibemMat(iP_x+l,iPxi+m) = integralEq14(iP_x,l,iPxi,m)
!     use resultVars, only: P => allpoints, B => boupoints
!                    
!     use Gquadrature, only: Gquad_n
!     implicit none
!     complex*16 :: integralEq16,desp
!     integer, intent(in) :: iP_x,l,iPxi,m
!     logical :: lejos
!     integer :: iGq
!     
!     integralEq16 = 0
!     lejos = .false.
!     
!     if (lejos) then
!     ! no usamos integración gaussiana
!            
!     else ! cerca
!     ! usamos integracion gaussiana
!     
!      ! En cada punto gaussiano del segmento en Xi en la dirección l
!      do iGq = 1,Gquad_n
!!        print*,"iGq=",iGq
!        ! dada la fuerza en XI en la dirección (m+1)
!        ! desplazamiento en la dirección (l+1)
!        !     receptor        emisor    
!        desp = P(iP_x)%GT_gq(iPxi,iGq,2-l,m+1)  
!!        print*,"desp=",desp
!        ! multiplicar por peso gaussiano de integración en segmento xi 
!        desp = desp * B(iPxi)%Gq_xXx_C(iGq)
!!        print*,"peso=",B(iPxi)%Gq_xXx_C(iGq)
!        ! acumular suma
!        integralEq16 = integralEq16 + desp
!      end do !xXx
!     end if !lejos
!!     print*,integralEq16 
!     end function integralEq16
!     
!     ! los deplazamientos
!     ! iMec     U        W
!     !  i       2        1
!     !  l       0        1
!     !necesito  x        z
!     function integralEq16mov(iP,iP_x,l,iPxi,m)
!        !                              receptor---,       ,--- fuente
!        !                                        _|___   _|___
!        !  ibemMat(iP_x+l,iPxi+m) = integralEq14(iP_x,l,iPxi,m)
!     use resultVars, only:P => allpoints, B => boupoints
!     use Gquadrature, only: Gquad_n
!     implicit none
!     complex*16 :: integralEq16mov,desp
!     integer, intent(in) :: iP,iP_x,l,iPxi,m
!     logical :: lejos
!     integer :: iGq
!     
!     integralEq16mov = 0
!     lejos = .false.
!     
!     if (lejos) then
!     ! no usamos integración gaussiana
!            
!     else ! cerca
!     ! usamos integracion gaussiana
!     
!      ! En cada punto gaussiano del segmento en Xi en la dirección l
!      do iGq = 1,Gquad_n
!        ! dada la fuerza en XI en la dirección (m+1)
!        ! desplazamiento en la dirección (l+1)
!        !     receptor        emisor    
!        desp = P(iP_x)%GT_gq_mov(iPxi,iGq,2-l,m+1,iP) 
!        ! multiplicar por peso gaussiano de integración en segmento xi 
!        desp = desp * B(iPxi)%Gq_xXx_C(iGq)
!        ! acumular suma
!        integralEq16mov = integralEq16mov + desp
!      end do !xXx
!     end if !lejos
!     
!     end function integralEq16mov
!     
!     ! los deplazamientos
!     ! iMec     U        W
!     !  i       2        1
!     !  l       0        1
!     !necesito  x        z
      subroutine W_to_t(W,iP,coords,outpf)
      use waveNumVars, only : NFREC,DFREC
      use glovars
      use waveVars, only : dt,t0,maxtime!,Uo
      use ploteo10pesos
      use wavelets
      implicit none
      integer ,intent(in) :: iP,outpf
      complex*16, dimension(Nfrec+1, imecMax), intent(in) :: W
      complex*16, dimension(2*NFREC) :: S
      real, dimension(2), intent(in) :: coords  
!     character(LEN=10),intent(in) :: tt
      character(LEN=3), dimension(4)   :: nombre 
      character(LEN=100) :: titleN
      character(LEN=9)   :: logflag
      integer :: i,iMec,t,n_maxtime
      real :: x_i,z_i,factor
      factor = sqrt(2.*NFREC)
      
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 'Tz-'
      nombre(4)= 'Tx-'
      
      x_i=coords(1)
      z_i=coords(2)
      
      do iMec = 1,imecMax
!     print*,sum(W(:,imec))
      !  (0) crepa en f
      !           1 2 3 4 ... Nfrec Nfrec+1
      ! tenemos:  0 1 2 3 ... N/2-1  N/2
      ! hacemos:  0 1 2 3 ... N/2-1 -N/2 ... -3 -2 -1
      S(      1:nfrec+1  ) =       W(1:nfrec+1:+1,  iMec) * & 
      exp(cmplx(0.0,- 2*pi * (/((t-1)*dfrec,t=1,nfrec+1)/) * t0,8))/(2*pi)!* Uo(1:nfrec+1)
      
      S(1) = W(1,iMec) * & 
      exp(cmplx(0.0,- 2*pi * 0.01*dfrec * t0,8))/(2*pi)!*0.01!* Uo(1)
      
!     S(1) = 0
      
      S(nfrec+2:nfrec*2) = conjg(W(nfrec:2:-1,iMec)) * & 
      exp(cmplx(0.0,- 2*pi * (/((-t)*dfrec,t=nfrec-1,1,-1)/) * t0,8))/(2*pi)!* Uo(nfrec+2:nfrec*2)
!     print*,"here 1"
      
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
        
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
!        print*,titleN
         call plotXYcomp(S,real(DFREC,4),2*nfrec,titleN, & 
         'frec[hz] ','amplitude',1200,800)
         
         if (Verbose .ge. 2) then
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'fL_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         logflag = 'logx     '
!        print*,"here1.5"
         end if
         call plotSpectrum(S,real(DFREC,4),2*nfrec,nfrec, & 
                         titleN,'frec[hz] ','amplitude',logflag,1200,800)
      !  (1) pasar al tiempo
         S = S*factor
         call fork(2*nfrec,S,+1,verbose,outpf)
         S = S/factor
      !  (2) remover efecto de la frecuencia imaginaria
         S = S * exp(DFREC/periodicdamper * Dt*((/(i,i=0,2*nfrec-1)/)))
         !tiempo maximo para graficar
         if(maxtime .lt. dt) maxtime = 5*dt
         if(maxtime .gt. 1/(dfrec)) maxtime = 1/(real(dfrec,4))
         n_maxtime = int(maxtime/real(dt,4))
      if (verbose .ge. 2) print*,"maxtime = ",maxtime," segs :: @",dt, & 
                                  " : ",n_maxtime," puntos"
      !  (3) plot the damn thing
      if (verbose .ge. 2) then
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
         do i = 1,n_maxtime
          write(3211,'(ES14.5E2,2x,ES14.5E2)') real(S(i)),aimag(S(i))
         end do
         close (3211) 
      end if
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'S_',nombre(iMec),iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S(1:n_maxtime),real(dt,4),n_maxtime,titleN, & 
         'time[sec]','amplitude',1200,800)
      
      end do !imec
      end subroutine W_to_t
      subroutine Hollywood(outpf)
      use DISLIN
      use resultVars
      use wavelets ! FORK(LX,CX,SIGNI,verbose,outpf)
      use waveVars, only : Dt,t0,maxtime!,Uo
      use waveNumVars ! FK,NFREC,DFREC
!     use refSolMatrixVars, only : COME
      use gloVars
      use MeshVars
      use soilVars, only : Z,N
      use GeometryVars, only : Xcoord,nxi
      implicit none
      integer ,intent(in) :: outpf
      real       ::  X(nxXxmpts)
      real       ::           Y(nz)
      complex*16 :: Sm(nxXxmpts,nz,imecMax+plusOne,2*NFREC)
      integer  :: iP,iMec
      integer  :: ix,iz,iT,i,Iframe,Fframe,mecMax,t, n_maxtime
      real     :: factor,ColorRangeMaximumScale,tlabel
      integer*4 :: lentitle
      character(LEN=3), dimension(5) :: nombre
      character(LEN=60) :: textoTimeStamp
      character(LEN=15) :: auxTxt
      CHARACTER(len=400) :: path
      character(LEN=1) :: imdone
      character(LEN=60) :: CTIT

      real :: maxY, minY, maxX, minX, p
      real, dimension(41)   :: ZLVRAY
      real                  :: maV,miV,Vstep,xstep,zstep
      real, parameter :: sharp = 14.
      real, dimension(imecMax+plusOne,2) :: colorBounds
      
      
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 'Tz-' 
      nombre(4)= 'Tx-'
      if (plotStress) then
      nombre(5)= 'mod'
      else
      nombre(3)= 'mod'
      end if
      mecMax = imecMax ! 2 o 4
      
      !tiempo maximo para graficar
         if(maxtime .lt. dt) maxtime = 5*dt
         if(maxtime .gt. 1/(dfrec)) maxtime = 1/(real(dfrec,4))
         n_maxtime = int(maxtime/dt,4)
         print*,"maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," puntos"
      
      if (verbose >= 1) Write(outpf,'(a)') "Will make a movie..."
         call system("mkdir video")
         CALL chdir("video")
         CALL getcwd(path)
         write(outpf,'(a,a)') "at ",TRIM(path)
         call system("rm *.png")
         call system("rm *.txt")
      
      factor = sqrt(real(2.*NFREC))
      
      allocate(Sismogram(size(X),size(Y),2,2*NFREC))
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
      do iP= mPtini,mPtfin ! cada profundidad de Z 
      do ix = 1,NxXxMpts   ! todos los pixeles en ese Z
        Y(iz) = allpoints(IP)%center(2) !;print*,iP,'z=',Y(iz)
        
        ! los espectros de cada coordenada ya están ordenaditos
        do iMec = 1,mecMax ! a 2 o 4
        ! (1) crepa
        Sm(ix,iz,iMec,1:nfrec+1) = allpoints(iP)%WmovieSiblings(ix,1:nfrec+1:+1,iMec)* & 
      exp(cmplx(0.0,- 2*pi*(/((t-1)*dfrec,t=1,nfrec+1)/) * t0,8))/(2*pi) !* Uo(1:nfrec+1)
      
        Sm(ix,iz,iMec,1      ) = allpoints(iP)%WmovieSiblings(ix,1         ,iMec) * & 
      exp(cmplx(0.0,- 2*pi * 0.01*dfrec * t0,8))/(2*pi)!*0.01 !* Uo(1)
      
        Sm(ix,iz,iMec,nfrec+2:nfrec*2) = conjg(allpoints(iP)%WmovieSiblings(ix,nfrec:2:-1,iMec))* & 
      exp(cmplx(0.0,- 2*pi*(/((-t)*dfrec,t=nfrec-1,1,-1)/) * t0,8))/(2*pi) !* Uo(nfrec+2:nfrec*2)
        
        ! (2) pasar al tiempo: fork
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) * factor
          call fork(2*nfrec,Sm(ix,iz,iMec,:),+1,verbose,outpf)
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) / factor
        
        ! (2.3) remover efecto de la frecuencia imaginaria
          Sm(ix,iz,iMec,:) = Sm(ix,iz,iMec,:) * & 
                            exp(DFREC/periodicdamper * Dt*((/(i,i=0,2*nfrec-1)/)))
          
          
        end do !iMec
!       ! el módulo de los desplazameintos
        if (plotModule) then
        Sm(ix,iz,imecMax+plusOne,:) = sqrt(Sm(ix,iz,1,:)**2 + Sm(ix,iz,2,:)**2)
        
        Sm(ix,iz,imecMax+plusOne,:) = log(1. + exp(sharp) * abs(Sm(ix,iz,imecMax+plusOne,:)))/ &
                        log(exp(sharp) + 1.)
        end if
      end do !ix
      iz = iz + 1
      end do !iP
        if (plotModule) then
        Sm(:,:,imecMax+plusOne,:) = Sm(:,:,imecMax+plusOne,:) / maxval(real(Sm(:,:,imecMax+plusOne,:)))
        end if
      !color table boundaries
      ColorRangeMaximumScale = 0.1
      
  123 do i=1,imecMax
       maV = maxVal(real(Sm(:,:,i,:),4))
       miV = minVal(real(Sm(:,:,i,:),4))
       maV = max(maV,abs(miV))
       miV = - maV
       
       if (verbose .ge. 2) then
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
       end if ! verbose
       
       colorBounds(i,1) = maV * ColorRangeMaximumScale
       colorBounds(i,2) = miV * ColorRangeMaximumScale
       
       if (verbose .ge. 2) then
          write(outpf,'(a,a,a,a,E12.4,a,E12.4,a)') "colorbounds:", & 
          "(",nombre(i),") [",colorBounds(i,2)," ; ",colorBounds(i,1),"]"
       end if
      end do !i
      ! y para el módulo
      if (plotModule) then
        colorBounds(imecMax+plusOne,1) = 1.
        colorBounds(imecMax+plusOne,2) = 0.
      end if

      minX=0.;maxX=0.;minY=0.;maxY=0.
      ! plotting boundaries according to geometry of topography
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
      Fframe = n_maxtime !size(Sm(1,1,1,:))
      
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
      do iMec=1,imecMax+plusOne !cada dirección y el módulo
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
      CALL PAGE (int(3100,4),int(2400,4))
      call imgfmt('RGB')
      call winsiz(int(1100,4),int(800,4)) !1200,800
      CALL SCRMOD('REVERS') !fondo blanco
      CALL DISINI()
!     CALL COMPLX ! sets a complex font
      CALL BMPFNT ('SIMPLEX')
!     CALL HWFONT()
           !the position of an axis system.
      CALL axspos (int(300,4) ,int(2200,4)) ! Lower left corner
      call axslen (int(2000,4), int(2000,4)) !size of the axis system.
      
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
      CALL SETVLT ('SPEC') !color table
!     call proj3d('ORTHO')
!     call view3d(real((maxX-minX)/2,4),real((maxY-minY)/2,4),real(1.1*maV,4),'ABS')
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4),& 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
      !CALL SHDMOD ('CELL', 'CONTUR')
      !CALL SHDMOD ('MIDDLE', 'COLOR')
      
      
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
      tlabel = (it)*real(dt,4)
!     print*,tlabel
      write(CTIT,'(a,F9.5,a)') 't=',tlabel,' seg'
      lentitle = NLMESS(CTIT)
      CALL MESSAG(CTIT,int((2900-lentitle-100),4),int(300,4))
      
      CALL HEIGHT(int(150,4))
      call errmod ("all", "off")
      CALL DISFIN 
      end do !iMec
!     stop "killed video"
      end do !iT
      
      !now make the video with the frames
      do iMec=1,imecmax+plusone
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
      end if!
      if (verbose .ge. 2) then
        write(6,'(a)') 'Delete *.txt files? [Y]'
        read(5,*)imdone
        if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
          write(path,'(a)') 'rm *.txt'
          call system(trim(path))
        end if
      end if
      end subroutine Hollywood

