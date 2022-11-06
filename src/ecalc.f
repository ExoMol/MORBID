      SUBROUTINE ECALC ( IWRKST , NREC2 ,
     1                   SKPSYM , MXJVP1 , ISO )
C*********************************************************************
C
C       PROGRAM FOR CALCULATING THE ROTATION-VIBRATION ENERGIES OF
C       A TRIATOMIC MOLECULE AFTER THE MORBID SCHEME OF PER JENSEN
C
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER ( LWORK = 6000000 )
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,12)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 , NFIL6 ,
     5       NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(55) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,12)
C
      LOGICAL SYMM,POTSYM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW,JDIMST(4),ISTART(4),
     1        LINDEX(4)
      LOGICAL SKPSYM ( 2 , 2 , MXJVP1 , ISOMAX )
      CHARACTER*126 ELINE
      CHARACTER*1 CABU
      CHARACTER*6 ROUTIN
      CHARACTER*3 INFOXX
      CHARACTER*4 SYMSYM(4),ASYSYM(4)
      INTEGER IWRKST,IWRKEN,INTWRK,
     1       NRECS,NRSYM0,NRASY0,NRSYMJ,NRASYJ,
     2       ISSYM0,ISASY0,ISSYMJ,ISASYJ,IENRGY,
     3       IF1ST,IF2ST,IV0ST,IRRST,IGFST,IPHIST,IDERST,
     4       IJACST,IRHSST,IXTXST,IXTYST,ICSSST,IASSST,IIA1ST,
     5       IIA2ST,IRTIRR,LENREC,       IVALST,IPHIL,IPHIR,
     6       IDERR,IDERL,ISTOST,NPHI,LENPHI,LENVAL,IIW1ST,IIW3ST,
     7       IOVEST,IVMAST,IHMAST,NOVA,NOVB,IAMAST,IEVMST,
     8       IUMAST,IEWKST,I,NFCT,LENIAR,LENINT,IRTIST,
     9       ICOMST,IEVIST,IDOMST,IDIM,NORECS,
     1       ITMAST,IESTST,ISMAST,IEWRST,INDBST,NAMDIM,
     1       IHTRST,IWMAST,V2VALU
      REAL*8 RECS
C
      COMMON /WORKCO/ WORK( LWORK )
C
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /DIMEN/  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
      DATA SYMSYM/'A1  ','B2  ','B1  ','A2  '/
      DATA ASYSYM/'A''  ','    ','A"  ','    '/
      CALL CLOCKV ( OLDVEC , OLDTIM , 1 , 2 )
C
C*********************************************************************
C
C       SEQUENCE FOR CALCULATING THE RHO DEPENDENT FUNCTIONS
C       NECESSARY FOR THE MORBID CALCULATION.
C
C*********************************************************************
C
      INTWRK=LWORK-IWRKST+1
C
C*********************************************************************
C
C THE RHO DEPENDENT FUNCTIONS ARE STORED ON DISK AS FOLLOWS
C
C NFIL1: FUNCTIONS NECESSARY TO DO A J=0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL2: ADDITIONAL FUNCTIONS NECESSARY TO DO A J=0 CALCULATION
C        FOR AN ABC MOLECULE.
C NFIL3: FUNCTIONS NECESSARY TO DO A J>0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL4: ADDITIONAL FUNCTIONS NECESSARY TO DO A J>0 CALCULATION
C        FOR AN ABC MOLECULE.
C
C*********************************************************************
C
      NFSYM0=31
      NFASY0=23
      NFSYMJ=52
      NFASYJ=40
      INTSYM=58
      INTASY=49
C
C*********************************************************************
C
C THE NF VARIABLES CONTAIN THE NUMBER OF FUNCTIONS WHOSE VALUES
C ARE TO BE STORED IN NFIL1, NFIL2, NFIL3, AND NFIL4, RESPECTIVELY.
C
C WE NOW DIVIDE UP THE WORK ARRAY TO HOLD THE VARIOUS VARIABLES
C
C*********************************************************************
C
      ISSYM0=IWRKST
      ISASY0=ISSYM0+NFSYM0
      ISSYMJ=ISASY0+NFASY0
      ISASYJ=ISSYMJ+NFSYMJ
      IWRKEN=ISASYJ+NFASYJ
      IF (IWRKEN .GT. LWORK) GOTO 1001
C
C*********************************************************************
C
C RHOFCT CALCULATES THE APPROPRIATE RHO DEPENDENT FUNCTIONS AND
C STORES THEM ON UNITS NFIL1, NFIL2, NFIL3, AND NFIL4.
C
C*********************************************************************
C
      CALL RHOFCT ( WORK(ISSYM0) , WORK(ISASY0) , WORK(ISSYMJ) ,
     1             WORK(ISASYJ) )
      CALL PRTIME( 'RHOFCT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C WE NOW RESERVE SPACE IN THE WORK ARRAY FOR THE PURE BENDING
C ENERGIES AND FOR THE FUNCTIONS NECESSARY FOR THE NUMEROV-COOLEY
C NUMERICAL INTEGRATION. THESE FUNCTIONS ARE: F1, F2, V0+UBEND,
C IRR0, AND G. FOR DEFINITIONS SEE P. JENSEN, COMP. PHYS. REP. 1,
C 1-55 (1983).
C
C*********************************************************************
C
      IENRGY=IWRKST
      IF1ST=IENRGY+V2MXP1*JMAXP1
      IF2ST=IF1ST+NSTINT
      IV0ST=IF2ST+NSTINT
      IRRST=IV0ST+NSTINT
      IGFST=IRRST+NSTINT
      IWRKEN=IGFST+NSTINT
      IF (IWRKEN .GE. LWORK) GOTO 1002
C
C*********************************************************************
C
C NUMFCT CALCULATES THE VALUES OF THESE FUNCTIONS
C
C*********************************************************************
C
      CALL NUMFCT ( WORK(IF1ST) , WORK(IF2ST) , WORK(IV0ST) ,
     1             WORK(IRRST) , WORK(IGFST) )
      CALL PRTIME( 'NUMFCT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C APART FROM THE SPACE ALREADY RESERVED IN THE WORK ARRAY
C WE NOW RESERVE SPACE FOR ONE BENDING WAVEFUNCTION (THE ONE
C THAT IS CURRENTLY BEING CALCULATED), FOR ITS DERIVATIVE
C AND FOR VARIOUS ARRAYS USED IN CARRYING OUT THE SERIES
C SOLUTION AROUND RHO=0.
C
C*********************************************************************
C
      IRTIST=IWRKEN
      IPHIST=IRTIST+NSTINT
      IDERST=IPHIST+NSTINT
      NSEPP2=NSERP+2
      NSEQP1=NSERQ+1
      IJACST=IDERST+NSTINT
      IRHSST=IJACST+NSERIN*NSEPP2
      IXTXST=IRHSST+NSERIN
      IXTYST=IXTXST+NSEPP2*NSEPP2
      ICSSST=IXTYST+NSEPP2
      IASSST=ICSSST+NSEPP2
      IIA1ST=IASSST+NSEQP1
      LENIAR=NSEPP2
      IIA2ST=IIA1ST+LENIAR
      IWRKEN=IIA2ST+LENIAR
      IF (IWRKEN .GT. LWORK) GOTO 1003
C
C*********************************************************************
C
C WAVFUN CALCULATES THE K=0 BENDING FUNCTIONS AND STORES THEM ON
C UNIT NFIL7.
C
C*********************************************************************
C
      CALL WAVFUN ( WORK(IENRGY) , WORK(IF1ST)  , WORK(IF2ST)  ,
     1             WORK(IV0ST)   , WORK(IRRST)  , WORK(IGFST)  ,
     2             WORK(IRTIST)  , WORK(IPHIST) , WORK(IDERST) ,
     3             WORK(IJACST)  , WORK(IRHSST) , WORK(IXTXST) ,
     4             WORK(IXTYST)  , WORK(ICSSST) , WORK(IASSST) ,
     5             WORK(IIA1ST)  , WORK(IIA2ST) ,
     6             LENIAR )
      CALL PRTIME( 'WAVFUN' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C WE NOW SET UP FOR CALCULATING THE INTEGRALS OVER THE BENDING
C COORDINATE NECESSARY FOR DOING A J=0 CALCULATION.
C WE SET UP SPACE FOR
C                       - THE INTEGRALS.
C                       - AS MANY VALUES OF THE RHO DEPENDENT
C                         THAT ONE NEEDS AT ANY GIVEN TIME.
C                       - A LEFT AND A RIGHT WAVEFUNCTION AND
C                         THEIR DERIVATIVES.
C
C*********************************************************************
C
      MBASP1=MBASIS+1
      NOBAS=MBASP1*(MBASP1+1)/2
      IESTST=IPHIST
      IVALST=IESTST+NOBAS
      IVSYST=IVALST+NSTINT*NFSYM0
      IF (NFSYM0.LT.NFASY0 .AND. .NOT.SYMM)
     1                     IVSYST=IVALST+NSTINT*NFASY0
      IPHIL=IVSYST+NSTINT*NFSYMJ
      IF (NFSYMJ .LT. NFASYJ .AND. .NOT. SYMM)
     1                       IPHIL=IVSYST+NSTINT*NFASYJ
      IDERL=IPHIL+NSTINT
      IPHIR=IDERL+NSTINT
      IDERR=IPHIR+NSTINT
      IWRKEN=IDERR+NSTINT
      LENVAL=IVSYST-IVALST
      LENVSY=IPHIL-IVSYST
      IF (IWRKEN .GT. LWORK) GOTO 1004
C
C*********************************************************************
C
C FIRST WE DO THE DELTA K=0 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTZER ( WORK(IGFST) , WORK(IVALST) , WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4                      LENVAL , LENVSY )
      CALL PRTIME( 'INTZER' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C THEN WE DO THE DELTA K=1 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTONE ( WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4             LENVSY )
      CALL PRTIME( 'INTONE' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C FINALLY WE DO THE DELTA K=2 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTTWO ( WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4             LENVSY )
      CALL PRTIME( 'INTTWO' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C SET UP PARAMETERS FOR THE CALCULATION OF MORSE OSCILLATOR MATRIX
C ELEMENTS
C
C*********************************************************************
C
      LENIAR=NOBAS
      MDIM=V2MXP1*(LSTYPA(1)+LSTYPA(2))
      IIW1ST=IVALST
      IIW3ST=IIW1ST+LENIAR
      IOVEST=IIW3ST+LENIAR
      LENIW=LENIAR
      IWRKEN=IOVEST+2*MBASP1*MBASP1
      IF (IWRKEN .GE. LWORK) GOTO 1005
C
C*********************************************************************
C
      CALL MOPARM ( WORK(IIW1ST) , WORK(IIW3ST) , WORK(IOVEST) )
      CALL PRTIME( 'MOPARM' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
      ITMAST=IWRKEN
      ISMAST=ITMAST+NOBAS*NOBAS
      IVMAST=ISMAST+NOBAS*NOBAS
      IUMAST=IVMAST+NOBAS*NOBAS
      IEWRST=IUMAST+NOBAS*2
      IHTRST=IEWRST+NOBAS
      ICOMST=IHTRST+NOBAS
      INDBST=ICOMST+LENIAR
      IDOMST=INDBST+LENIAR
      IWRKEN=IDOMST+LENIAR
      IF (IWRKEN .GE. LWORK) GOTO 1006
C
C*********************************************************************
C
      CALL STRTCH ( WORK(IIW1ST) , WORK(IIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IUMAST) , WORK(IEWRST) ,
     1             WORK(IHTRST) , WORK(ICOMST) , WORK(INDBST) ,
     1             WORK(IDOMST) )
      CALL PRTIME( 'STRTCH' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
      NOPERA=38
      NSTMAT=NOPERA*NOBAS*NOBAS
      IWMAST=IUMAST
      IWRKEN=IWMAST+NOBAS*NOBAS
      NOPERA=38
      NSTMAT=NOPERA*NOBAS*NOBAS
      KSTMAT=IIW1ST
      IIW1ST=IIW1ST+NSTMAT
      IIW3ST=IIW3ST+NSTMAT
      IOVEST=IOVEST+NSTMAT
      ITMAST=ITMAST+NSTMAT
      ISMAST=ISMAST+NSTMAT
      IVMAST=IVMAST+NSTMAT
      IWMAST=IWMAST+NSTMAT
      IWRKEN=IWRKEN+NSTMAT
      IF (IWRKEN .GE. LWORK) GOTO 1007
      DO 12 III=IIW1ST,IWRKEN-1
12    WORK(III)=WORK(III-NSTMAT)
C
C*********************************************************************
C
      CALL GENSTR( WORK(IIW1ST) , WORK(IIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IWMAST) , WORK(KSTMAT) ,
     1             NOPERA )
      CALL PRTIME( 'GENSTR' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C  START LOOP OVER J = 0 THROUGH JMAX
C
C*********************************************************************
C
      ISYMST=IIW1ST
      IASYST=ISYMST+INTSYM
      ISMAST=IASYST+INTASY
      IROTOT=(2*JMAX+1)*V2MXP1*(LSTYPA(1)+LSTYPA(2))
      ISTORE=ISMAST+NOBAS*NOBAS
      ISTOR2=ISTORE+IROTOT
      IHMAST=ISTOR2+IROTOT
C
      DO 200 JP1=1,JMAXP1
C
      J=JP1-1
      INDEXJ=0
      INDEXD=1
      DO 160 II=1,4
160   JDIMST(II)=0
C
      DO 180 ITAUP1=1,2
      ITAU=ITAUP1-1
C
      IF (MOD(J,2) .EQ. 0) THEN
            IROTE=J/2+1
            IROTO=J/2
      ELSE
            IROTE=(J+1)/2
            IROTO=IROTE
      ENDIF
      IF (MOD(J+ITAU,2) .NE. 0) IROTE=IROTE-1
C
C*********************************************************************
C
      IF (.NOT. SYMM) THEN
            IRODIM=(IROTE+IROTO)*V2MXP1*LSTYPA(1)
      ELSE
            IRODIM=(IROTE*LSTYPA(1)+IROTO*LSTYPA(2))*V2MXP1
      ENDIF
      IWRKEN=IHMAST+IRODIM*IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1010
C
C*********************************************************************
C
C  CLEAR ALL ELEMENTS IN THE HAMILTONIAN MATRIX-TO-BE
C
C*********************************************************************
C
      DO 100 IH=IHMAST,IWRKEN
100   WORK(IH)=0.0D+00
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (SKPSYM(ITAUP1,1,JP1,ISO)) GOTO 380
C
      IF (IROTE .EQ. 0) GOTO 320
      IEO=0
      NOFSET=0
      NDIMST=LSTYPA(1)
      ISTOFF=0
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , INTSYM , INTASY , ISTOFF ,
     1             WORK(KSTMAT) , NOPERA )
      NOFSRO=0
      NOFSCO=0
      NDIMRO=LSTYPA(1)
      NDIMCO=LSTYPA(1)
      ISTOFF=0
      JSTOFF=0
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
C
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
320   CONTINUE
      IF (IROTO .EQ. 0) GOTO 340
      IEO=1
      NOFSET=IROTE*V2MXP1*NDIMST
      IF (SYMM) THEN
            NDIMST=LSTYPA(2)
            ISTOFF=KSTYPA(1)
      ENDIF
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , INTSYM , INTASY , ISTOFF ,
     1             WORK(KSTMAT) , NOPERA )
      NOFSRO=NOFSET
      NOFSCO=NOFSET
      IF (SYMM) THEN
            NDIMRO=LSTYPA(2)
            NDIMCO=LSTYPA(2)
            ISTOFF=KSTYPA(1)
            JSTOFF=KSTYPA(1)
      ENDIF
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
C
C*********************************************************************
C
C  SET UP DELTA K = 1 MATRIX ELEMENTS CONNECTING THE TWO CORNERS
C
C*********************************************************************
C
340   CONTINUE
      IF (IROTO .EQ. 0 .OR. IROTE .EQ. 0) GOTO 360
      IEO=0
      NOFSRO=0
      NOFSCO=NOFSET
      IF (SYMM) THEN
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(2)
            ISTOFF=0
            JSTOFF=KSTYPA(1)
      ENDIF
      CALL DELK1 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
360   CONTINUE
C
C*********************************************************************
C
C  DIAGONALIZE THE HAMILTONIAN MATRIX BLOCK AND PRINT OUT THE
C  EIGENVALUES
C
C*********************************************************************
C
      IF (IRODIM .EQ. 0) GOTO 380
      IEIGST=IWRKEN
      IEWKST=IEIGST+IRODIM
      IQNUMS=IEWKST+IRODIM
      IWRKEN=IQNUMS+IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1011
      CALL DIAROT ( IRODIM , WORK(IHMAST) , WORK(IEIGST) ,
     1 WORK(IEWKST) , WORK(IQNUMS) ,
     1 J  ,         ITAU   , 1            )
C
      IF (JP1 .EQ. 1 .AND. ITAU .EQ. 0) THEN
            E00=WORK(IEIGST)
            IF (IPRINT .GT. 0) WRITE (NFIL6,6600) E00
      ENDIF
      DO 165 II=1,IRODIM
      WORK(ISTORE+INDEXJ)=WORK(IEIGST+II-1)
      WORK(ISTOR2+INDEXJ)=WORK(IQNUMS+II-1)
      WORK(II           )=WORK(IEIGST+II-1)
165   INDEXJ=INDEXJ+1
      JDIMST(INDEXD)=IRODIM
C
380   INDEXD=INDEXD+1
C
      DO 7000 II=IRODIM+1,NREC2
7000  WORK(II)=0.0D+00
      IREC=4*J+2*ITAU+1
      WRITE (NFIL16,REC=IREC) (WORK(II), II=1,NREC2)
C
C*********************************************************************
C
C
      IF ((.NOT. SYMM) .OR. SKPSYM(ITAUP1,2,JP1,ISO)) GOTO 460
C
C*********************************************************************
C
      IRODIM=(IROTE*LSTYPA(2)+IROTO*LSTYPA(1))*V2MXP1
      IWRKEN=IHMAST+IRODIM*IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1010
C
C*********************************************************************
C
      DO 120 IH=IHMAST,IWRKEN
120   WORK(IH)=0.0D+00
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (IROTE .EQ. 0) GOTO 400
      IEO=0
      NOFSET=0
      NDIMST=LSTYPA(2)
      ISTOFF=KSTYPA(1)
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , INTSYM , INTASY , ISTOFF ,
     1             WORK(KSTMAT) , NOPERA )
      NOFSRO=0
      NOFSCO=0
      NDIMRO=LSTYPA(2)
      NDIMCO=LSTYPA(2)
      ISTOFF=KSTYPA(1)
      JSTOFF=KSTYPA(1)
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
C
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
400   CONTINUE
      IF (IROTO .EQ. 0) GOTO 420
      IEO=1
      NOFSET=IROTE*V2MXP1*NDIMST
      NDIMST=LSTYPA(1)
      ISTOFF=0
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , INTSYM , INTASY , ISTOFF ,
     1             WORK(KSTMAT) , NOPERA )
      NOFSRO=NOFSET
      NOFSCO=NOFSET
      NDIMRO=LSTYPA(1)
      NDIMCO=LSTYPA(1)
            ISTOFF=0
            JSTOFF=0
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
C
C*********************************************************************
C
C  SET UP DELTA K = 1 MATRIX ELEMENTS CONNECTING THE TWO CORNERS
C
C*********************************************************************
C
420   CONTINUE
      IF (IROTE .EQ. 0 .OR. IROTO .EQ. 0) GOTO 440
      IEO=0
      NOFSRO=0
      NOFSCO=NOFSET
      NDIMRO=LSTYPA(2)
      NDIMCO=LSTYPA(1)
      ISTOFF=KSTYPA(1)
      JSTOFF=0
      CALL DELK1 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , INTSYM ,
     4             INTASY , ISTOFF , JSTOFF ,
     1             WORK(KSTMAT) , NOPERA )
440   CONTINUE
C
C*********************************************************************
C
      IF (IRODIM .EQ. 0) GOTO 460
      IEIGST=IWRKEN
      IEWKST=IEIGST+IRODIM
      IQNUMS=IEWKST+IRODIM
      IWRKEN=IQNUMS+IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1011
      CALL DIAROT ( IRODIM , WORK(IHMAST) , WORK(IEIGST) ,
     1 WORK(IEWKST) , WORK(IQNUMS) ,
     1 J  ,         ITAU   , 2            )
C
      DO 170 II=1,IRODIM
      WORK(ISTORE+INDEXJ)=WORK(IEIGST+II-1)
      WORK(ISTOR2+INDEXJ)=WORK(IQNUMS+II-1)
      WORK(II           )=WORK(IEIGST+II-1)
170   INDEXJ=INDEXJ+1
      JDIMST(INDEXD)=IRODIM
460   INDEXD=INDEXD+1
C
      DO 7010 II=IRODIM+1,NREC2
7010  WORK(II)=0.0D+00
      IREC=4*J+2*ITAU+2
      WRITE (NFIL16,REC=IREC) (WORK(II), II=1,NREC2)
C
C*********************************************************************
C
180   CONTINUE
C
      IF (IPRINT .LE. 0) GOTO 200
C
      IOFFS = ISTOR2 - ISTORE
C
      WRITE (NFIL6,6500) J
C
      ISTART(1)=ISTORE
      ISTART(2)=ISTART(1)+JDIMST(1)
      ISTART(3)=ISTART(2)+JDIMST(2)
      ISTART(4)=ISTART(3)+JDIMST(3)
C
      DO 182 II=1,4
182   LINDEX(II)=0
C
      II=0
      DO 184 JJ=1,4
      IF (JDIMST(JJ) .NE. 0) THEN
            II=II+1
            LINDEX(II)=JJ
      ENDIF
184   CONTINUE
      NCOL=II
C
      DO 290 II=1,126
290   ELINE(II:II)=' '
      DO 300 II=1,11
300   ELINE(II:II)='-'
      DO 620 JJ=1,NCOL
      IPLLO=12+(JJ-1)*28
      IPLHI=IPLLO+27
      DO 310 II=IPLLO,IPLHI
310   ELINE(II:II)='-'
620   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 330 II=1,126
330   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      DO 640 II=1,NCOL
      JJ=39+(II-1)*28
640   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 350 II=1,NCOL
      JJ=24+(II-1)*28
      KINDEX=LINDEX(II)
      IF (SYMM) THEN
            WRITE (ELINE(JJ:JJ+3),6400) SYMSYM(KINDEX)
      ELSE
            WRITE (ELINE(JJ:JJ+3),6400) ASYSYM(KINDEX)
      ENDIF
350   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 660 II=1,126
660   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      DO 680 II=1,NCOL
      JJ=39+(II-1)*28
680   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 661 II=1,126
661   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      DO 681 II=1,NCOL
      JJ=39+(II-1)*28
      WRITE (ELINE(JJ-27:JJ-1),6111)
681   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 390 II=1,11
390   ELINE(II:II)='-'
      DO 410 JJ=1,NCOL
      IPLLO=12+(JJ-1)*28
      IPLHI=IPLLO+27
      DO 700 II=IPLLO,IPLHI
700   ELINE(II:II)='-'
410   CONTINUE
      WRITE (NFIL6,6300) ELINE
C
      IDIMMX=MAX0(JDIMST(1),JDIMST(2),JDIMST(3),JDIMST(4))
C
      DO 412 II=1,126
412   ELINE(II:II)=' '
      DO 190 II=1,IDIMMX
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      KINDEX=II-1
      WRITE (ELINE(4:6),6000) KINDEX
C
      DO 188 JJ=1,NCOL
      KINDEX=LINDEX(JJ)
      IPLLO=12+(JJ-1)*28
      IPLHI=IPLLO+25
      IF (JDIMST(KINDEX) .GE. II) THEN
           EOUT = WORK(ISTART(KINDEX)+II-1) - E00
           XPACK= WORK(ISTART(KINDEX)+II-1 + IOFFS)
           XZ=XPACK/1.0D+05
           KVALU=INT(XZ)
           XPACK=XPACK - KVALU*1.0D+05
           XZ=XPACK/1.0D+03
           NVALU=INT(XZ)
           XPACK=XPACK - NVALU*1.0D+03
           XZ=XPACK/1.0D+01
           V2VALU=INT(XZ)
           XPACK=XPACK - V2VALU*1.0D+01
           IABU=NINT(XPACK)
           IF (IABU .EQ. 1) THEN
                CABU='A'
           ELSE
                CABU='B'
           ENDIF
           WRITE (ELINE(IPLLO:IPLHI),6100) KVALU,V2VALU,NVALU,
     1     CABU,EOUT
      ELSE
           WRITE (ELINE(IPLLO:IPLHI),6200)
      ENDIF
188   ELINE(IPLHI+2:IPLHI+2)=':'
190   WRITE (NFIL6,6300) ELINE
      DO 720 II=1,11
720   ELINE(II:II)='-'
      DO 740 JJ=1,NCOL
      IPLLO=12+(JJ-1)*28
      IPLHI=IPLLO+27
      DO 430 II=IPLLO,IPLHI
430   ELINE(II:II)='-'
740   CONTINUE
      WRITE (NFIL6,6300) ELINE
      CALL PRTIMJ( J , OLDTIM , OLDVEC )
200   CONTINUE
10    CONTINUE
6000  FORMAT(I3)
6100  FORMAT(3I4,1X,A1,1X,F11.5)
6111  FORMAT('  KA  V2  NS         E     ')
6200  FORMAT(27X)
6300  FORMAT(1H ,1X,A126)
6400  FORMAT(A4)
6500  FORMAT(1H0,12('*'),'  J = ',I3,' ENERGIES  ',12('*'),/)
6600  FORMAT(1H0,12('*'),'  ZERO POINT ENERGY  ',12('*'),//,
     1      12X,'E0 = ',F12.5//)
      RETURN
1001  WRITE (NFIL6,2001)
2001  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' INITIAL FUNCTION CALCULATION')
      STOP
1002  WRITE (NFIL6,2002)
2002  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATION OF F1, F2 ETC.')
      STOP
1003  WRITE (NFIL6,2003)
2003  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' NUMEROV-COOLEY INTEGRATION')
      STOP
1004  WRITE (NFIL6,2004)
2004  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' BENDING INTEGRAL CALCULATION (J=0)')
      STOP
1005  WRITE (NFIL6,2005)
2005  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' SETTING UP MORSE OSCILLATOR PARAMETERS')
      STOP
1006  WRITE (NFIL6,2006)
2006  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' DIAGONALIZING THE STRETCHING MATRIX')
      STOP
1007  WRITE (NFIL6,2007)
2007  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATING STRETCHING MATRIX ELEMENTS')
      STOP
1010  WRITE (NFIL6,2010) J
2010  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' ROTATIONAL CALCULATION AT J = ',I3)
      STOP
1011  WRITE (NFIL6,2011) J
2011  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' DIAGONALIZATON AT J = ',I3)
      END
C
C
      SUBROUTINE INTZER ( GRHO   , FCTVAL , ROTVAL , PHIL   ,
     1                   DERL   , PHIR   , DERR   ,
     3                            LENVAL , LENVSY )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,12)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5       NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(55) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,12)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER I,L,LENINT,LENVAL,LENREC,NPHI,INDRC ,
     1       NORECS , ILEN ,ISHFT
      REAL*8 GRHO(NSTINT),               FCTVAL(LENVAL),PHIL(NSTINT),
     1      DERL(NSTINT),PHIR(NSTINT),DERR(NSTINT),ROTVAL(LENVSY)
      REAL*8  SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,SUM08,SUM09,
     1 SUM10,SUM11,SUM12,SUM13,SUM14,SUM15,SUM16,SUM17,SUM18,SUM19,
     2 SUM20,SUM21,SUM22,SUM23,SUM24,SUM25,SUM26,SUM27,SUM28,SUM29,
     3 SUM30,SUM31
      REAL*8 FACODD,FACEVE
      INTEGER INDEX1,INDEX2,
     1       IOFFS,IOFFS0,IEND,JOFFS,KOFFSO,KOFFSE,NFCT
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /DIMEN/  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9001)
      FACODD=2.0D+00*HSTEP/3.0D+00
      FACEVE=4.0D+00*HSTEP/3.0D+00
C
C READ IN FUNCTION VALUES FROM DISK
C
      REWIND NFIL1
      REWIND NFIL3
      DO 1099 I=1,NSTINT,2
      IOFFS1=(I-1)*NFSYM0
      IOFFS2=(I-1)*NFSYMJ
      READ (NFIL1) (FCTVAL(IOFFS1+L), L=1,NFSYM0)
      READ (NFIL1) (FCTVAL(IOFFS1+L+NFSYM0), L=1,NFSYM0)
      READ (NFIL3) (ROTVAL(IOFFS2+L), L=1,NFSYMJ)
      READ (NFIL3) (ROTVAL(IOFFS2+L+NFSYMJ), L=1,NFSYMJ)
1099  CONTINUE
C
C START LOOP OVER K FOR CALCULATING DELTAK=0 MATRIX ELEMENTS
C
      DO 1410 KP1=1,JMAXP1
      KSQUAR=(KP1-1)*(KP1-1)
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 1400 INDEX1=1,V2MXP1
C
      IREC1=(KP1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 1200 INDEX2=1,V2MXP1
C
      IREC2=(KP1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM01=0.0D+00
      SUM02=0.0D+00
      SUM03=0.0D+00
      SUM04=0.0D+00
      SUM05=0.0D+00
      SUM06=0.0D+00
      SUM07=0.0D+00
      SUM08=0.0D+00
      SUM09=0.0D+00
      SUM10=0.0D+00
      SUM11=0.0D+00
      SUM12=0.0D+00
      SUM13=0.0D+00
      SUM14=0.0D+00
      SUM15=0.0D+00
      SUM16=0.0D+00
      SUM17=0.0D+00
      SUM18=0.0D+00
      SUM19=0.0D+00
      SUM20=0.0D+00
      SUM21=0.0D+00
      SUM22=0.0D+00
      SUM23=0.0D+00
      SUM24=0.0D+00
      SUM25=0.0D+00
      SUM26=0.0D+00
      SUM27=0.0D+00
      SUM28=0.0D+00
      SUM29=0.0D+00
      SUM30=0.0D+00
      SUM31=0.0D+00
      SUM32=0.0D+00
      SUM33=0.0D+00
      SUM34=0.0D+00
      SUM35=0.0D+00
      SUM36=0.0D+00
      SUM37=0.0D+00
      SUM38=0.0D+00
      SUM39=0.0D+00
      SUM40=0.0D+00
      SUM41=0.0D+00
      SUM42=0.0D+00
      SUM43=0.0D+00
      SUM44=0.0D+00
      SUM45=0.0D+00
      SUM46=0.0D+00
      SUM47=0.0D+00
      SUM48=0.0D+00
      SUM49=0.0D+00
      SUM50=0.0D+00
      SUM51=0.0D+00
      SUM52=0.0D+00
      SUM53=0.0D+00
      SUM54=0.0D+00
      SUM55=0.0D+00
      SUM56=0.0D+00
      SUM57=0.0D+00
      SUM58=0.0D+00
C
      DO 1100 I=1,NSTINT,2
      IOFFS1=(I-1)*NFSYM0
      IOFFS2=(I-1)*NFSYMJ
      SUM01=SUM01+FACODD*FCTVAL(IOFFS1+ 1)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 1+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM02=SUM02+FACODD*FCTVAL(IOFFS1+ 2)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 2+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM03=SUM03+FACODD*FCTVAL(IOFFS1+ 3)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 3+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM04=SUM04+FACODD*FCTVAL(IOFFS1+ 4)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 4+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM05=SUM05+FACODD*FCTVAL(IOFFS1+ 5)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 5+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM06=SUM06+FACODD*FCTVAL(IOFFS1+ 6)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 6+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM07=SUM07+FACODD*FCTVAL(IOFFS1+ 7)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 7+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM08=SUM08+FACODD*FCTVAL(IOFFS1+ 8)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 8+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM09=SUM09+FACODD*FCTVAL(IOFFS1+ 9)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 9+NFSYM0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM10=SUM10+FACODD*PHIL(I  )*FCTVAL(IOFFS1+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+10+NFSYM0)*PHIR(I+1)
      SUM11=SUM11+FACODD*PHIL(I  )*FCTVAL(IOFFS1+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+11+NFSYM0)*PHIR(I+1)
      SUM12=SUM12+FACODD*PHIL(I  )*FCTVAL(IOFFS1+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+12+NFSYM0)*PHIR(I+1)
      SUM13=SUM13+FACODD*PHIL(I  )*FCTVAL(IOFFS1+13)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+13+NFSYM0)*PHIR(I+1)
      SUM14=SUM14+FACODD*PHIL(I  )*FCTVAL(IOFFS1+14)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+14+NFSYM0)*PHIR(I+1)
      SUM15=SUM15+FACODD*PHIL(I  )*FCTVAL(IOFFS1+15)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+15+NFSYM0)*PHIR(I+1)
      SUM16=SUM16+FACODD*PHIL(I  )*FCTVAL(IOFFS1+16)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+16+NFSYM0)*PHIR(I+1)
      SUM17=SUM17+FACODD*PHIL(I  )*FCTVAL(IOFFS1+17)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+17+NFSYM0)*PHIR(I+1)
      SUM18=SUM18+FACODD*PHIL(I  )*FCTVAL(IOFFS1+18)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+18+NFSYM0)*PHIR(I+1)
      SUM19=SUM19+FACODD*PHIL(I  )*FCTVAL(IOFFS1+19)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+19+NFSYM0)*PHIR(I+1)
      SUM20=SUM20+FACODD*PHIL(I  )*FCTVAL(IOFFS1+20)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+20+NFSYM0)*PHIR(I+1)
      SUM21=SUM21+FACODD*PHIL(I  )*FCTVAL(IOFFS1+21)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+21+NFSYM0)*PHIR(I+1)
      SUM22=SUM22+FACODD*PHIL(I  )*FCTVAL(IOFFS1+22)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+22+NFSYM0)*PHIR(I+1)
      SUM23=SUM23+FACODD*PHIL(I  )*FCTVAL(IOFFS1+23)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+23+NFSYM0)*DERR(I+1)
      SUM24=SUM24+FACODD*PHIL(I  )*FCTVAL(IOFFS1+24)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+24+NFSYM0)*DERR(I+1)
      SUM25=SUM25+FACODD*PHIL(I  )*FCTVAL(IOFFS1+25)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+25+NFSYM0)*DERR(I+1)
      SUM26=SUM26+FACODD*PHIL(I  )*FCTVAL(IOFFS1+26)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+26+NFSYM0)*DERR(I+1)
      SUM27=SUM27+FACODD*PHIL(I  )*FCTVAL(IOFFS1+27)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+27+NFSYM0)*DERR(I+1)
      SUM28=SUM28+FACODD*PHIL(I  )*FCTVAL(IOFFS1+28)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+28+NFSYM0)*DERR(I+1)
      SUM29=SUM29+FACODD*PHIL(I  )*FCTVAL(IOFFS1+29)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+29+NFSYM0)*DERR(I+1)
      SUM30=SUM30+FACODD*PHIL(I  )*FCTVAL(IOFFS1+30)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+30+NFSYM0)*DERR(I+1)
      SUM31=SUM31+FACODD*PHIL(I  )*FCTVAL(IOFFS1+31)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+31+NFSYM0)*DERR(I+1)
      SUM32=SUM32+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 1)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 1+NFSYMJ)*PHIR(I+1)
      SUM33=SUM33+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 2)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 2+NFSYMJ)*PHIR(I+1)
      SUM34=SUM34+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 3)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 3+NFSYMJ)*PHIR(I+1)
      SUM35=SUM35+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 4)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 4+NFSYMJ)*PHIR(I+1)
      SUM36=SUM36+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 5)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 5+NFSYMJ)*PHIR(I+1)
      SUM37=SUM37+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 6)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 6+NFSYMJ)*PHIR(I+1)
      SUM38=SUM38+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 7)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 7+NFSYMJ)*PHIR(I+1)
      SUM39=SUM39+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 8)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 8+NFSYMJ)*PHIR(I+1)
      SUM40=SUM40+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 9)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 9+NFSYMJ)*PHIR(I+1)
      SUM41=SUM41+FACODD*PHIL(I  )*ROTVAL(IOFFS2+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+10+NFSYMJ)*PHIR(I+1)
      SUM42=SUM42+FACODD*PHIL(I  )*ROTVAL(IOFFS2+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+11+NFSYMJ)*PHIR(I+1)
      SUM43=SUM43+FACODD*PHIL(I  )*ROTVAL(IOFFS2+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+12+NFSYMJ)*PHIR(I+1)
      SUM44=SUM44+FACODD*PHIL(I  )*ROTVAL(IOFFS2+13)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+13+NFSYMJ)*PHIR(I+1)
      SUM45=SUM45+FACODD*PHIL(I  )*ROTVAL(IOFFS2+14)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+14+NFSYMJ)*PHIR(I+1)
      SUM46=SUM46+FACODD*PHIL(I  )*ROTVAL(IOFFS2+15)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+15+NFSYMJ)*PHIR(I+1)
      SUM47=SUM47+FACODD*PHIL(I  )*ROTVAL(IOFFS2+16)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+16+NFSYMJ)*PHIR(I+1)
      SUM48=SUM48+FACODD*PHIL(I  )*ROTVAL(IOFFS2+17)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+17+NFSYMJ)*PHIR(I+1)
      SUM49=SUM49+FACODD*PHIL(I  )*ROTVAL(IOFFS2+18)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+18+NFSYMJ)*PHIR(I+1)
      SUM50=SUM50+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+19)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+19+NFSYMJ)
     2          *PHIR(I+1)
      SUM51=SUM51+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+20)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+20+NFSYMJ)
     2          *PHIR(I+1)
      SUM52=SUM52+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+21)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+21+NFSYMJ)
     2          *PHIR(I+1)
      SUM53=SUM53+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+22)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+22+NFSYMJ)
     2          *PHIR(I+1)
      SUM54=SUM54+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+23)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+23+NFSYMJ)
     2          *PHIR(I+1)
      SUM55=SUM55+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+24)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+24+NFSYMJ)
     2          *PHIR(I+1)
      SUM56=SUM56+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+25)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+25+NFSYMJ)
     2          *PHIR(I+1)
      SUM57=SUM57+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+26)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+26+NFSYMJ)
     2          *PHIR(I+1)
      SUM58=SUM58+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+27)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+27+NFSYMJ)
     2          *PHIR(I+1)
1100  CONTINUE
      WRITE (NFIL9) SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,
     1              SUM08,SUM09,SUM10,SUM11,SUM12,SUM13,SUM14,
     2              SUM15,SUM16,SUM17,SUM18,SUM19,SUM20,SUM21,
     3              SUM22,SUM23,SUM24,SUM25,SUM26,SUM27,SUM28,
     4              SUM29,SUM30,SUM31,SUM32,SUM33,SUM34,SUM35,
     1              SUM36,SUM37,SUM38,SUM39,SUM40,SUM41,SUM42,
     5              SUM43,SUM44,SUM45,SUM46,SUM47,SUM48,SUM49,
     6              SUM50,SUM51,SUM52,SUM53,SUM54,SUM55,SUM56,
     7              SUM57,SUM58
1200  CONTINUE
1400  CONTINUE
1410  CONTINUE
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9002)
      IF (SYMM) RETURN
      IOFFS0=NFSYM0*V2MXP1*V2MXP1
C
C*********************************************************************
C
C THIS IS THE END OF THE CALCULATION FOR THE SYMMETRICAL MOLECULE.
C WE NOW CALCULATE INTEGRALS FOR THE UNSYMMETRICAL CASE.
C
C*********************************************************************
C
C  WE FIRST READ IN THE NECESSARY MATRIX ELEMENTS
C
      REWIND NFIL2
      REWIND NFIL4
      DO 2809 I=1,NSTINT,2
      IOFFS1=(I-1)*NFASY0
      IOFFS2=(I-1)*NFASYJ
      READ (NFIL2) (FCTVAL(IOFFS1+L), L=1,NFASY0)
      READ (NFIL2) (FCTVAL(IOFFS1+L+NFASY0), L=1,NFASY0)
      READ (NFIL4) (ROTVAL(IOFFS2+L), L=1,NFASYJ)
      READ (NFIL4) (ROTVAL(IOFFS2+L+NFASYJ), L=1,NFASYJ)
2809  CONTINUE

C
C START LOOP OVER K FOR CALCULATING THE DELTAK=0 MATRIX ELEMENTS
C
      DO 2810 KP1=1,JMAXP1
      KSQUAR=(KP1-1)*(KP1-1)
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 2800 INDEX1=1,V2MXP1
C
      IREC1=(KP1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 2600 INDEX2=1,V2MXP1
      IREC2=(KP1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM01=0.0D+00
      SUM02=0.0D+00
      SUM03=0.0D+00
      SUM04=0.0D+00
      SUM05=0.0D+00
      SUM06=0.0D+00
      SUM07=0.0D+00
      SUM08=0.0D+00
      SUM09=0.0D+00
      SUM10=0.0D+00
      SUM11=0.0D+00
      SUM12=0.0D+00
      SUM13=0.0D+00
      SUM14=0.0D+00
      SUM15=0.0D+00
      SUM16=0.0D+00
      SUM17=0.0D+00
      SUM18=0.0D+00
      SUM19=0.0D+00
      SUM20=0.0D+00
      SUM21=0.0D+00
      SUM22=0.0D+00
      SUM23=0.0D+00
      SUM24=0.0D+00
      SUM25=0.0D+00
      SUM26=0.0D+00
      SUM27=0.0D+00
      SUM28=0.0D+00
      SUM29=0.0D+00
      SUM30=0.0D+00
      SUM31=0.0D+00
      SUM32=0.0D+00
      SUM33=0.0D+00
      SUM34=0.0D+00
      SUM35=0.0D+00
      SUM36=0.0D+00
      SUM37=0.0D+00
      SUM38=0.0D+00
      SUM39=0.0D+00
      SUM40=0.0D+00
      SUM41=0.0D+00
C
      DO 2500 I=1,NSTINT,2
      IOFFS1=(I-1)*NFASY0
      IOFFS2=(I-1)*NFASYJ
      SUM01=SUM01+FACODD*FCTVAL(IOFFS1+ 1)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 1+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM02=SUM02+FACODD*FCTVAL(IOFFS1+ 2)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 2+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM03=SUM03+FACODD*FCTVAL(IOFFS1+ 3)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 3+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM04=SUM04+FACODD*FCTVAL(IOFFS1+ 4)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 4+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM05=SUM05+FACODD*FCTVAL(IOFFS1+ 5)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 5+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM06=SUM06+FACODD*FCTVAL(IOFFS1+ 6)*(DERL(I  )*DERR(I  )
     1                          +GRHO(I  )*PHIL(I  )*PHIR(I  ))
     2          +FACEVE*FCTVAL(IOFFS1+ 6+NFASY0)*(DERL(I+1)*DERR(I+1)
     3                          +GRHO(I+1)*PHIL(I+1)*PHIR(I+1))
      SUM07=SUM07+FACODD*PHIL(I  )*FCTVAL(IOFFS1+ 7)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+ 7+NFASY0)*PHIR(I+1)
      SUM08=SUM08+FACODD*PHIL(I  )*FCTVAL(IOFFS1+ 8)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+ 8+NFASY0)*PHIR(I+1)
      SUM09=SUM09+FACODD*PHIL(I  )*FCTVAL(IOFFS1+ 9)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+ 9+NFASY0)*PHIR(I+1)
      SUM10=SUM10+FACODD*PHIL(I  )*FCTVAL(IOFFS1+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+10+NFASY0)*PHIR(I+1)
      SUM11=SUM11+FACODD*PHIL(I  )*FCTVAL(IOFFS1+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+11+NFASY0)*PHIR(I+1)
      SUM12=SUM12+FACODD*PHIL(I  )*FCTVAL(IOFFS1+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+12+NFASY0)*PHIR(I+1)
      SUM13=SUM13+FACODD*PHIL(I  )*FCTVAL(IOFFS1+13)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+13+NFASY0)*PHIR(I+1)
      SUM14=SUM14+FACODD*PHIL(I  )*FCTVAL(IOFFS1+14)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+14+NFASY0)*PHIR(I+1)
      SUM15=SUM15+FACODD*PHIL(I  )*FCTVAL(IOFFS1+15)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+15+NFASY0)*DERR(I+1)
      SUM16=SUM16+FACODD*PHIL(I  )*FCTVAL(IOFFS1+16)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+16+NFASY0)*DERR(I+1)
      SUM17=SUM17+FACODD*PHIL(I  )*FCTVAL(IOFFS1+17)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+17+NFASY0)*DERR(I+1)
      SUM18=SUM18+FACODD*PHIL(I  )*FCTVAL(IOFFS1+18)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+18+NFASY0)*DERR(I+1)
      SUM19=SUM19+FACODD*PHIL(I  )*FCTVAL(IOFFS1+19)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+19+NFASY0)*DERR(I+1)
      SUM20=SUM20+FACODD*PHIL(I  )*FCTVAL(IOFFS1+20)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+20+NFASY0)*DERR(I+1)
      SUM21=SUM21+FACODD*PHIL(I  )*FCTVAL(IOFFS1+21)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+21+NFASY0)*DERR(I+1)
      SUM22=SUM22+FACODD*PHIL(I  )*FCTVAL(IOFFS1+22)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+22+NFASY0)*DERR(I+1)
      SUM23=SUM23+FACODD*PHIL(I  )*FCTVAL(IOFFS1+23)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*FCTVAL(IOFFS1+23+NFASY0)*DERR(I+1)
      SUM24=SUM24+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 1)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 1+NFASYJ)*PHIR(I+1)
      SUM25=SUM25+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 2)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 2+NFASYJ)*PHIR(I+1)
      SUM26=SUM26+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 3)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 3+NFASYJ)*PHIR(I+1)
      SUM27=SUM27+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 4)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 4+NFASYJ)*PHIR(I+1)
      SUM28=SUM28+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 5)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 5+NFASYJ)*PHIR(I+1)
      SUM29=SUM29+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 6)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 6+NFASYJ)*PHIR(I+1)
      SUM30=SUM30+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 7)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 7+NFASYJ)*PHIR(I+1)
      SUM31=SUM31+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 8)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 8+NFASYJ)*PHIR(I+1)
      SUM32=SUM32+FACODD*PHIL(I  )*ROTVAL(IOFFS2+ 9)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+ 9+NFASYJ)*PHIR(I+1)
      SUM33=SUM33+FACODD*PHIL(I  )*ROTVAL(IOFFS2+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+10+NFASYJ)*PHIR(I+1)
      SUM34=SUM34+FACODD*PHIL(I  )*ROTVAL(IOFFS2+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+11+NFASYJ)*PHIR(I+1)
      SUM35=SUM35+FACODD*PHIL(I  )*ROTVAL(IOFFS2+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS2+12+NFASYJ)*PHIR(I+1)
      SUM36=SUM36+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+13)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+13+NFASYJ)
     1          *PHIR(I+1)
      SUM37=SUM37+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+14)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+14+NFASYJ)
     1          *PHIR(I+1)
      SUM38=SUM38+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+15)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+15+NFASYJ)
     1          *PHIR(I+1)
      SUM39=SUM39+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+16)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+16+NFASYJ)
     1          *PHIR(I+1)
      SUM40=SUM40+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+17)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+17+NFASYJ)
     1          *PHIR(I+1)
      SUM41=SUM41+FACODD*PHIL(I  )*KSQUAR*ROTVAL(IOFFS2+18)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*KSQUAR*ROTVAL(IOFFS2+18+NFASYJ)
     1          *PHIR(I+1)
2500  CONTINUE
      WRITE (NFIL10) SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,
     1              SUM08,SUM09,SUM10,SUM11,SUM12,SUM13,SUM14,
     2              SUM15,SUM16,SUM17,SUM18,SUM19,SUM20,SUM21,
     3              SUM22,SUM23,SUM24,SUM25,SUM26,SUM27,SUM28,
     4              SUM29,SUM30,SUM31,SUM32,SUM33,SUM34,SUM35,
     5              SUM36,SUM37,SUM38,SUM39,SUM40,SUM41
2600  CONTINUE
2800  CONTINUE
2810  CONTINUE
      IF (IPRINT .NE. 0) WRITE (NFIL6,9003)
      RETURN
9001  FORMAT(1H0,' MORBID.INT.INF   DK=0 INTEGRATION ROUTINE ENTERED'/)
9002  FORMAT(1H0,' MORBID.INT.INF   SYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=0'/)
9003  FORMAT(1H0,' MORBID.INT.INF   UNSYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=0'/)
      END
C
C
      SUBROUTINE INTONE ( ROTVAL , PHIL   ,
     1                   DERL   , PHIR   , DERR   ,
     3                   LENVSY )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,12)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5       NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(55) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,12)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      REAL*8 PHIL(NSTINT),
     1      DERL(NSTINT),PHIR(NSTINT),DERR(NSTINT),ROTVAL(LENVSY)
      REAL*8  SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,SUM08,SUM09,
     1 SUM10,SUM11,SUM12,SUM13,SUM14,SUM15,SUM16,SUM17,SUM18,SUM19,
     2 SUM20,SUM21,SUM22,SUM23,SUM24,SUM25,SUM26,SUM27,SUM28,SUM29,
     3 SUM30,SUM31
      REAL*8 FACODD,FACEVE
      INTEGER INDEX1,INDEX2,
     1       IOFFS,IOFFS0,IEND,JOFFS,KOFFSO,KOFFSE,NFCT
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /DIMEN/  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9001)
      FACODD=2.0D+00*HSTEP/3.0D+00
      FACEVE=4.0D+00*HSTEP/3.0D+00
C
C     READ NECESSARY MATRIX ELEMENTS
C
      REWIND NFIL3
      DO 1099 I=1,NSTINT,2
      IOFFS=(I-1)*NFSYMJ
      READ (NFIL3) (ROTVAL(IOFFS+L), L=1,NFSYMJ)
      READ (NFIL3) (ROTVAL(IOFFS+L+NFSYMJ), L=1,NFSYMJ)
1099  CONTINUE
C
C START LOOP OVER K FOR CALCULATING DELTAK=1 MATRIX ELEMENTS
C
      DO 1410 K1P1=1,JMAXP1
      DO 1405 KDEL=-1,1,2
      K2P1=K1P1+KDEL
      IF (K2P1 .LT. 1 .OR. K2P1 .GT. JMAXP1) GOTO 1405
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 1400 INDEX1=1,V2MXP1
C
      IREC1=(K1P1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 1200 INDEX2=1,V2MXP1
C
      IREC2=(K2P1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM01=0.0D+00
      SUM02=0.0D+00
      SUM03=0.0D+00
      SUM04=0.0D+00
      SUM05=0.0D+00
      SUM06=0.0D+00
      SUM07=0.0D+00
      SUM08=0.0D+00
      SUM09=0.0D+00
      SUM10=0.0D+00
      SUM11=0.0D+00
      SUM12=0.0D+00
      SUM13=0.0D+00
      SUM14=0.0D+00
      SUM15=0.0D+00
      SUM16=0.0D+00
      SUM17=0.0D+00
      SUM18=0.0D+00
      SUM19=0.0D+00
      SUM20=0.0D+00
      SUM21=0.0D+00
      SUM22=0.0D+00
      SUM23=0.0D+00
      SUM24=0.0D+00
      SUM25=0.0D+00
C
      DO 1100 I=1,NSTINT,2
      IOFFS=(I-1)*NFSYMJ
      SUM01=SUM01+FACODD*PHIL(I  )*ROTVAL(IOFFS+28)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+28+NFSYMJ)*PHIR(I+1)
      SUM02=SUM02+FACODD*PHIL(I  )*ROTVAL(IOFFS+29)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+29+NFSYMJ)*PHIR(I+1)
      SUM03=SUM03+FACODD*PHIL(I  )*ROTVAL(IOFFS+30)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+30+NFSYMJ)*PHIR(I+1)
      SUM04=SUM04+FACODD*PHIL(I  )*ROTVAL(IOFFS+31)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+31+NFSYMJ)*PHIR(I+1)
      SUM05=SUM05+FACODD*PHIL(I  )*ROTVAL(IOFFS+32)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+32+NFSYMJ)*PHIR(I+1)
      SUM06=SUM06+FACODD*PHIL(I  )*ROTVAL(IOFFS+33)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+33+NFSYMJ)*PHIR(I+1)
      SUM07=SUM07+FACODD*PHIL(I  )*ROTVAL(IOFFS+34)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+34+NFSYMJ)*PHIR(I+1)
      SUM08=SUM08+FACODD*PHIL(I  )*ROTVAL(IOFFS+35)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+35+NFSYMJ)*PHIR(I+1)
      SUM09=SUM09+FACODD*PHIL(I  )*ROTVAL(IOFFS+36)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+36+NFSYMJ)*DERR(I+1)
      SUM10=SUM10+FACODD*PHIL(I )*ROTVAL(IOFFS+37)*DERR(I )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+37+NFSYMJ)*DERR(I+1)
      SUM11=SUM11+FACODD*PHIL(I  )*ROTVAL(IOFFS+38)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+38+NFSYMJ)*DERR(I+1)
      SUM12=SUM12+FACODD*PHIL(I  )*ROTVAL(IOFFS+39)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+39+NFSYMJ)*DERR(I+1)
      SUM13=SUM13+FACODD*PHIL(I  )*ROTVAL(IOFFS+40)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+40+NFSYMJ)*DERR(I+1)
      SUM14=SUM14+FACODD*PHIL(I  )*ROTVAL(IOFFS+41)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+41+NFSYMJ)*DERR(I+1)
      SUM15=SUM15+FACODD*PHIL(I  )*ROTVAL(IOFFS+42)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+42+NFSYMJ)*DERR(I+1)
      SUM16=SUM16+FACODD*PHIL(I  )*ROTVAL(IOFFS+43)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+43+NFSYMJ)*DERR(I+1)
      SUM17=SUM17+FACODD*PHIL(I  )*ROTVAL(IOFFS+44)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+44+NFSYMJ)*PHIR(I+1)
      SUM18=SUM18+FACODD*PHIL(I  )*ROTVAL(IOFFS+45)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+45+NFSYMJ)*PHIR(I+1)
      SUM19=SUM19+FACODD*PHIL(I  )*ROTVAL(IOFFS+46)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+46+NFSYMJ)*PHIR(I+1)
      SUM20=SUM20+FACODD*PHIL(I  )*ROTVAL(IOFFS+47)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+47+NFSYMJ)*PHIR(I+1)
      SUM21=SUM21+FACODD*PHIL(I  )*ROTVAL(IOFFS+48)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+48+NFSYMJ)*PHIR(I+1)
      SUM22=SUM22+FACODD*PHIL(I  )*ROTVAL(IOFFS+49)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+49+NFSYMJ)*PHIR(I+1)
      SUM23=SUM23+FACODD*PHIL(I  )*ROTVAL(IOFFS+50)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+50+NFSYMJ)*PHIR(I+1)
      SUM24=SUM24+FACODD*PHIL(I  )*ROTVAL(IOFFS+51)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+51+NFSYMJ)*PHIR(I+1)
      SUM25=SUM25+FACODD*PHIL(I  )*ROTVAL(IOFFS+52)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+52+NFSYMJ)*PHIR(I+1)
1100  CONTINUE
      WRITE (NFIL11) SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,
     1              SUM08,SUM09,SUM10,SUM11,SUM12,SUM13,SUM14,
     2              SUM15,SUM16,SUM17,SUM18,SUM19,SUM20,SUM21,
     3              SUM22,SUM23,SUM24,SUM25
1200  CONTINUE
1400  CONTINUE
1405  CONTINUE
1410  CONTINUE
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9002)
      IF (SYMM) RETURN
      IOFFS0=NFSYM0*V2MXP1*V2MXP1
C
C*********************************************************************
C
C THIS THE THE END OF THE CALCULATION FOR THE SYMMETRICAL MOLECULE.
C WE NOW CALCULATE INTEGRALS FOR THE UNSYMMETRICAL CASE.
C
C*********************************************************************
C
C     READ NECESSARY MATRIX ELEMENTS
C
      REWIND NFIL4
      DO 2499 I=1,NSTINT,2
      IOFFS=(I-1)*NFASYJ
      READ (NFIL4) (ROTVAL(IOFFS+L), L=1,NFASYJ)
      READ (NFIL4) (ROTVAL(IOFFS+L+NFASYJ), L=1,NFASYJ)
2499  CONTINUE
C
C START LOOP OVER K FOR CALCULATING THE DELTA K=1 MATRIX ELEMENTS
C
      DO 2810 K1P1=1,JMAXP1
      DO 2805 KDEL=-1,1,2
      K2P1=K1P1+KDEL
      IF (K2P1 .LT. 1 .OR. K2P1 .GT. JMAXP1) GOTO 2805
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 2800 INDEX1=1,V2MXP1
C
      IREC1=(K1P1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 2600 INDEX2=1,V2MXP1
      IREC2=(K2P1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM01=0.0D+00
      SUM02=0.0D+00
      SUM03=0.0D+00
      SUM04=0.0D+00
      SUM05=0.0D+00
      SUM06=0.0D+00
      SUM07=0.0D+00
      SUM08=0.0D+00
      SUM09=0.0D+00
      SUM10=0.0D+00
      SUM11=0.0D+00
      SUM12=0.0D+00
      SUM13=0.0D+00
      SUM14=0.0D+00
      SUM15=0.0D+00
      SUM16=0.0D+00
      SUM17=0.0D+00
      SUM18=0.0D+00
      SUM19=0.0D+00
      SUM20=0.0D+00
      SUM21=0.0D+00
      SUM22=0.0D+00
C
      DO 2500 I=1,NSTINT,2
      IOFFS=(I-1)*NFASYJ
      SUM01=SUM01+FACODD*PHIL(I  )*ROTVAL(IOFFS+19)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+19+NFASYJ)*PHIR(I+1)
      SUM02=SUM02+FACODD*PHIL(I  )*ROTVAL(IOFFS+20)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+20+NFASYJ)*PHIR(I+1)
      SUM03=SUM03+FACODD*PHIL(I  )*ROTVAL(IOFFS+21)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+21+NFASYJ)*PHIR(I+1)
      SUM04=SUM04+FACODD*PHIL(I  )*ROTVAL(IOFFS+22)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+22+NFASYJ)*PHIR(I+1)
      SUM05=SUM05+FACODD*PHIL(I  )*ROTVAL(IOFFS+23)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+23+NFASYJ)*PHIR(I+1)
      SUM06=SUM06+FACODD*PHIL(I  )*ROTVAL(IOFFS+24)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+24+NFASYJ)*PHIR(I+1)
      SUM07=SUM07+FACODD*PHIL(I  )*ROTVAL(IOFFS+25)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+25+NFASYJ)*PHIR(I+1)
      SUM08=SUM08+FACODD*PHIL(I  )*ROTVAL(IOFFS+26)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+26+NFASYJ)*DERR(I+1)
      SUM09=SUM09+FACODD*PHIL(I  )*ROTVAL(IOFFS+27)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+27+NFASYJ)*DERR(I+1)
      SUM10=SUM10+FACODD*PHIL(I  )*ROTVAL(IOFFS+28)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+28+NFASYJ)*DERR(I+1)
      SUM11=SUM11+FACODD*PHIL(I  )*ROTVAL(IOFFS+29)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+29+NFASYJ)*DERR(I+1)
      SUM12=SUM12+FACODD*PHIL(I  )*ROTVAL(IOFFS+30)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+30+NFASYJ)*DERR(I+1)
      SUM13=SUM13+FACODD*PHIL(I  )*ROTVAL(IOFFS+31)*DERR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+31+NFASYJ)*DERR(I+1)
      SUM14=SUM14+FACODD*PHIL(I  )*ROTVAL(IOFFS+32)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+32+NFASYJ)*PHIR(I+1)
      SUM15=SUM15+FACODD*PHIL(I  )*ROTVAL(IOFFS+33)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+33+NFASYJ)*PHIR(I+1)
      SUM16=SUM16+FACODD*PHIL(I  )*ROTVAL(IOFFS+34)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+34+NFASYJ)*PHIR(I+1)
      SUM17=SUM17+FACODD*PHIL(I  )*ROTVAL(IOFFS+35)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+35+NFASYJ)*PHIR(I+1)
      SUM18=SUM18+FACODD*PHIL(I  )*ROTVAL(IOFFS+36)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+36+NFASYJ)*PHIR(I+1)
      SUM19=SUM19+FACODD*PHIL(I  )*ROTVAL(IOFFS+37)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+37+NFASYJ)*PHIR(I+1)
      SUM20=SUM20+FACODD*PHIL(I  )*ROTVAL(IOFFS+38)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+38+NFASYJ)*PHIR(I+1)
      SUM21=SUM21+FACODD*PHIL(I  )*ROTVAL(IOFFS+39)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+39+NFASYJ)*PHIR(I+1)
      SUM22=SUM22+FACODD*PHIL(I  )*ROTVAL(IOFFS+40)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+40+NFASYJ)*PHIR(I+1)
2500  CONTINUE
      WRITE (NFIL12) SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,
     1              SUM08,SUM09,SUM10,SUM11,SUM12,SUM13,SUM14,
     2              SUM15,SUM16,SUM17,SUM18,SUM19,SUM20,SUM21,
     3              SUM22
2600  CONTINUE
2800  CONTINUE
2805  CONTINUE
2810  CONTINUE
      IF (IPRINT .NE. 0) WRITE (NFIL6,9003)
      RETURN
9001  FORMAT(1H0,' MORBID.INT.INF   DK=1 INTEGRATION ROUTINE ENTERED'/)
9002  FORMAT(1H0,' MORBID.INT.INF   SYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=1'/)
9003  FORMAT(1H0,' MORBID.INT.INF   UNSYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=1'/)
      END
C
C
      SUBROUTINE INTTWO ( ROTVAL , PHIL   ,
     1                   DERL   , PHIR   , DERR   ,
     3                   LENVSY )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,12)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5       NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(55) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,12)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      REAL*8 PHIL(NSTINT),
     1      DERL(NSTINT),PHIR(NSTINT),DERR(NSTINT),ROTVAL(LENVSY)
      REAL*8  SUM01,SUM02,SUM03,SUM04,SUM05,SUM06,SUM07,SUM08,SUM09,
     1 SUM10,SUM11,SUM12,SUM13,SUM14,SUM15,SUM16,SUM17,SUM18,SUM19,
     2 SUM20,SUM21,SUM22,SUM23,SUM24,SUM25,SUM26,SUM27,SUM28,SUM29,
     3 SUM30,SUM31
      REAL*8 FACODD,FACEVE
      INTEGER INDEX1,INDEX2,
     1       IOFFS,IOFFS0,IEND,JOFFS,KOFFSO,KOFFSE,NFCT
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /DIMEN/  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9001)
      FACODD=2.0D+00*HSTEP/3.0D+00
      FACEVE=4.0D+00*HSTEP/3.0D+00
C
C     READ NECESSARY MATRIX ELEMENTS
C
      REWIND NFIL3
      DO 1099 I=1,NSTINT,2
      IOFFS=(I-1)*NFSYMJ
      READ (NFIL3) (ROTVAL(IOFFS+L), L=1,NFSYMJ)
      READ (NFIL3) (ROTVAL(IOFFS+L+NFSYMJ), L=1,NFSYMJ)
1099  CONTINUE
C
C START LOOP OVER K FOR CALCULATING DELTAK=2 MATRIX ELEMENTS
C
      DO 1410 K1P1=1,JMAXP1
      DO 1405 KDEL=-2,2,4
      K2P1=K1P1+KDEL
      IF (K2P1 .LT. 1 .OR. K2P1 .GT. JMAXP1) GOTO 1405
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 1400 INDEX1=1,V2MXP1
C
      IREC1=(K1P1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 1200 INDEX2=1,V2MXP1
C
      IREC2=(K2P1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM32=0.0D+00
      SUM33=0.0D+00
      SUM34=0.0D+00
      SUM35=0.0D+00
      SUM36=0.0D+00
      SUM37=0.0D+00
      SUM38=0.0D+00
      SUM39=0.0D+00
      SUM40=0.0D+00
      SUM41=0.0D+00
      SUM42=0.0D+00
      SUM43=0.0D+00
      SUM44=0.0D+00
      SUM45=0.0D+00
      SUM46=0.0D+00
      SUM47=0.0D+00
      SUM48=0.0D+00
      SUM49=0.0D+00
C
      DO 1100 I=1,NSTINT,2
      IOFFS=(I-1)*NFSYMJ
      SUM32=SUM32+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 1)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 1+NFSYMJ)*PHIR(I+1)
      SUM33=SUM33+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 2)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 2+NFSYMJ)*PHIR(I+1)
      SUM34=SUM34+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 3)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 3+NFSYMJ)*PHIR(I+1)
      SUM35=SUM35+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 4)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 4+NFSYMJ)*PHIR(I+1)
      SUM36=SUM36+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 5)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 5+NFSYMJ)*PHIR(I+1)
      SUM37=SUM37+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 6)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 6+NFSYMJ)*PHIR(I+1)
      SUM38=SUM38+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 7)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 7+NFSYMJ)*PHIR(I+1)
      SUM39=SUM39+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 8)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 8+NFSYMJ)*PHIR(I+1)
      SUM40=SUM40+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 9)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 9+NFSYMJ)*PHIR(I+1)
      SUM41=SUM41+FACODD*PHIL(I  )*ROTVAL(IOFFS+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+10+NFSYMJ)*PHIR(I+1)
      SUM42=SUM42+FACODD*PHIL(I  )*ROTVAL(IOFFS+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+11+NFSYMJ)*PHIR(I+1)
      SUM43=SUM43+FACODD*PHIL(I  )*ROTVAL(IOFFS+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+12+NFSYMJ)*PHIR(I+1)
      SUM44=SUM44+FACODD*PHIL(I  )*ROTVAL(IOFFS+13)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+13+NFSYMJ)*PHIR(I+1)
      SUM45=SUM45+FACODD*PHIL(I  )*ROTVAL(IOFFS+14)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+14+NFSYMJ)*PHIR(I+1)
      SUM46=SUM46+FACODD*PHIL(I  )*ROTVAL(IOFFS+15)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+15+NFSYMJ)*PHIR(I+1)
      SUM47=SUM47+FACODD*PHIL(I  )*ROTVAL(IOFFS+16)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+16+NFSYMJ)*PHIR(I+1)
      SUM48=SUM48+FACODD*PHIL(I  )*ROTVAL(IOFFS+17)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+17+NFSYMJ)*PHIR(I+1)
      SUM49=SUM49+FACODD*PHIL(I  )*ROTVAL(IOFFS+18)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+18+NFSYMJ)*PHIR(I+1)
1100  CONTINUE
      WRITE (NFIL13)       SUM32,SUM33,SUM34,SUM35,
     1              SUM36,SUM37,SUM38,SUM39,SUM40,SUM41,SUM42,
     5              SUM43,SUM44,SUM45,SUM46,SUM47,SUM48,SUM49
1200  CONTINUE
1400  CONTINUE
1405  CONTINUE
1410  CONTINUE
C
      IF (IPRINT .NE. 0) WRITE (NFIL6,9002)
      IF (SYMM) RETURN
      IOFFS0=NFSYM0*V2MXP1*V2MXP1
C
C*********************************************************************
C
C THIS THE THE END OF THE CALCULATION FOR THE SYMMETRICAL MOLECULE.
C WE NOW CALCULATE INTEGRALS FOR THE UNSYMMETRICAL CASE.
C
C*********************************************************************
C
C     READ NECESSARY MATRIX ELEMENTS
C
      REWIND NFIL4
      DO 2499 I=1,NSTINT,2
      IOFFS=(I-1)*NFASYJ
      READ (NFIL4) (ROTVAL(IOFFS+L), L=1,NFASYJ)
      READ (NFIL4) (ROTVAL(IOFFS+L+NFASYJ), L=1,NFASYJ)
2499  CONTINUE
C
C START LOOP OVER K FOR CALCULATING THE DELTAK=2 MATRIX ELEMENTS
C
      DO 2810 K1P1=1,JMAXP1
      DO 2805 KDEL=-2,2,4
      K2P1=K1P1+KDEL
      IF (K2P1 .LT. 1 .OR. K2P1 .GT. JMAXP1) GOTO 2805
C
C START LOOP OVER V2(LEFT) FOR THE CALCULATION OF MATRIX ELEMENTS
C
      DO 2800 INDEX1=1,V2MXP1
C
      IREC1=(K1P1-1)*V2MXP1+INDEX1
      READ (NFIL7,REC=IREC1) (PHIL(I),I=1,NSTINT),
     1                       (DERL(I),I=1,NSTINT)
C
      DO 2600 INDEX2=1,V2MXP1
      IREC2=(K2P1-1)*V2MXP1+INDEX2
      READ (NFIL7,REC=IREC2) (PHIR(I),I=1,NSTINT),
     1                       (DERR(I),I=1,NSTINT)
C
      SUM24=0.0D+00
      SUM25=0.0D+00
      SUM26=0.0D+00
      SUM27=0.0D+00
      SUM28=0.0D+00
      SUM29=0.0D+00
      SUM30=0.0D+00
      SUM31=0.0D+00
      SUM32=0.0D+00
      SUM33=0.0D+00
      SUM34=0.0D+00
      SUM35=0.0D+00
C
      DO 2500 I=1,NSTINT,2
      IOFFS=(I-1)*NFASYJ
      SUM24=SUM24+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 1)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 1+NFASYJ)*PHIR(I+1)
      SUM25=SUM25+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 2)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 2+NFASYJ)*PHIR(I+1)
      SUM26=SUM26+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 3)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 3+NFASYJ)*PHIR(I+1)
      SUM27=SUM27+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 4)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 4+NFASYJ)*PHIR(I+1)
      SUM28=SUM28+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 5)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 5+NFASYJ)*PHIR(I+1)
      SUM29=SUM29+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 6)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 6+NFASYJ)*PHIR(I+1)
      SUM30=SUM30+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 7)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 7+NFASYJ)*PHIR(I+1)
      SUM31=SUM31+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 8)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 8+NFASYJ)*PHIR(I+1)
      SUM32=SUM32+FACODD*PHIL(I  )*ROTVAL(IOFFS+ 9)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+ 9+NFASYJ)*PHIR(I+1)
      SUM33=SUM33+FACODD*PHIL(I  )*ROTVAL(IOFFS+10)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+10+NFASYJ)*PHIR(I+1)
      SUM34=SUM34+FACODD*PHIL(I  )*ROTVAL(IOFFS+11)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+11+NFASYJ)*PHIR(I+1)
      SUM35=SUM35+FACODD*PHIL(I  )*ROTVAL(IOFFS+12)*PHIR(I  )
     1          +FACEVE*PHIL(I+1)*ROTVAL(IOFFS+12+NFASYJ)*PHIR(I+1)
2500  CONTINUE
      WRITE (NFIL14) SUM24,SUM25,SUM26,SUM27,SUM28,
     4              SUM29,SUM30,SUM31,SUM32,SUM33,SUM34,SUM35
2600  CONTINUE
2800  CONTINUE
2805  CONTINUE
2810  CONTINUE
      IF (IPRINT .NE. 0) WRITE (NFIL6,9003)
      RETURN
9001  FORMAT(1H0,' MORBID.INT.INF   DK=2 INTEGRATION ROUTINE ENTERED'/)
9002  FORMAT(1H0,' MORBID.INT.INF   SYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=2'/)
9003  FORMAT(1H0,' MORBID.INT.INF   UNSYMMETRIC MOLECULE INTEGRALS CA',
     1            'LCULATED FOR DK=2'/)
      END

