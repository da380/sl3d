
C$**********************************************************************
      SUBROUTINE ROTMX2(NMAX,L,THETA,D,ID1,ID2)
c
c
c     This subroutine calculates the 'Wigner d-functions'
c     d_{nm}(theta) that arise from the matrix elements of 
c     irreducible unitary representations of the rotation 
c     group SO(3).
c
c     Inputs:
c
c            NMAX  = maximum value for the column index
c                    the input array D should have column dimension
c                    at least 2*NMAX+1
c
c           L      = the angular degree of the function 
c
c           THETA  = the rotation angle in radians
c
c           ID1    = column dimension of D
c
c           ID2    = row dimension of D: must be at least 2L+1
c
c      Outputs
c             D    = an array of of the d-functions of dimension 
c                   (2*NMAX+1)-by-(2*L+1) such that the element d_{nm}
c                   is  equal to D(n+NMAX+1,m+L+1).
c
c      The d_{nm} functions are such that d_{nm}(0) is the identity 
c      matrix, and are the same as the generalized Legendre polynomials
c      occuring in Phinney and Burridge. Because of this, the generalized
c      spherical harmonic functions Y^{N}_{lm} are given by
c
c           Y^{N}_{lm}(theta,phi)  = d_{nm}(theta)exp(i*m*phi)
c
c
c     Such generalized spherical harmonics are normalized such that
c     their L^{2} norm over the unit sphere is equal to (4*pi)/(2L+1)
c    
c
c     Fully normalized generalized spherical harmonics such as those 
c     occurring in Dahlen and Tromp can be formed from those above 
c     by multiplying them by sqrt((2L+1)/(4*pi))
c
c
c==================================================================%
c          This routine was written by John Woodhouse.             %         
c==================================================================%

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      save
      DOUBLE PRECISION D,THETA
      DIMENSION D(ID1,ID2)
      DATA BIG,SMALL,DLBIG,DLSML/1.D35,1.D-35,35.D0,-35.D0/
      DATA PI/3.141592653589793238462643383279502884197D0/
      DFLOAT(N)=N
      TH=THETA
      IF(TH.GT.PI) TH = PI
      IF(TH.LT.0.D0) TH = 0.D0
!      IF((TH.GT.PI).OR.(TH.LT.0.D0))  STOP 'ILLEGAL ARG IN ROTMX2'
      IF(L.NE.0) GOTO 350
      D(1+NMAX,L+1)=1.D0
      RETURN
350   ISUP=1
      IF(TH.LE.PI/2.D0) GOTO 310
      TH=PI-TH
      ISUP=-1
310   NM=2*L+1
      NMP1=NM+1
      LP1=L+1
      LM1=L-1
      LP2=L+2
      NROW=2*NMAX+1
      NMAXP1=NMAX+1
      LMN=L-NMAX
      IF(TH.NE.0.D0) GOTO 320
      DO 330 IM1CT=1,NROW
      IM1=IM1CT+LMN
      DO 330 IM2=LP1,NM
      D(IM1CT,IM2)=0.D0
      IF(IM1.EQ.IM2) D(IM1CT,IM2)=1.D0
330   CONTINUE
      GOTO 400
320   CONTINUE
C
C     ZERO L.H.S. OF MATRIX
C
      DO 340 IM1=1,NROW
      DO 340 IM2=1,LP1
340   D(IM1,IM2)=0.D0
C
C        SET UP PARAMETERS
C
      SHTH=DSIN(0.5D0*TH)
      CHTH=DCOS(0.5D0*TH)
      STH=2.D0*SHTH*CHTH
      CTH=2.D0*CHTH*CHTH-1.D0
      DLOGF=DLOG10(CHTH/SHTH)
      DLOGS=DLOG10(SHTH)
C
C       ITERATE FROM LAST COLUMN USING 1. AS STARTING VALUE
C
      DO 10 IM1CT=1,NROW
      IM1=IM1CT+LMN
      M1=IM1-LP1
      RM1=M1
      NM2=MIN0(IM1-1,NM-IM1)
      D(IM1CT,NM)=1.D0
      IF(NM2.EQ.0) GOTO 10
      DO 20 NIT=1,NM2
      M2=L-NIT
      IM2=M2+LP1
      IF(M2.NE.LM1) GOTO 70
      T1=0.D0
      GOTO 30
70    T1=-DSQRT(DFLOAT((IM2+1)*(L-M2-1)))*D(IM1CT,IM2+2)
30    D(IM1CT,IM2)=T1-(2.D0/STH)*(CTH*DFLOAT(M2+1)-RM1)
     1    *D(IM1CT,IM2+1)
      D(IM1CT,IM2)=D(IM1CT,IM2)/DSQRT(DFLOAT(IM2*(L-M2)))
      TEMP=D(IM1CT,IM2)
      RMOD=DABS(TEMP)
      IF(RMOD.LT.BIG) GOTO 20
      IF(NIT.EQ.NM2) GOTO 20
      D(IM1CT,NIT+1)=DLBIG
      D(IM1CT,IM2)=D(IM1CT,IM2)/BIG
      D(IM1CT,IM2+1)=D(IM1CT,IM2+1)/BIG
20    CONTINUE
10    CONTINUE
C
C        SET UP NORMALIZATION FOR RIGHTMOST COLUMN
C
      T1=DFLOAT(2*L)*DLOGS
      IF(LMN.EQ.0) GOTO 720
      DO 710 I=1,LMN
      M1=I-L
      T1=DLOGF+0.5D0*DLOG10(DFLOAT(LP1-M1)/DFLOAT(L+M1))+T1
710   CONTINUE
720   D(1,1)=T1
      IF(NROW.EQ.1) GOTO 730
      DO 110 IM1CT=2,NROW
      M1=IM1CT-NMAXP1
110   D(IM1CT,1)=DLOGF+0.5D0*DLOG10(DFLOAT(L-M1+1)/DFLOAT(L+M1))
     1     +D(IM1CT-1,1)
730   SGN=-1.D0
      IF((LMN/2)*2.NE.LMN) SGN=1.D0
C
C       RENORMALIZE ROWS
C
      DO 120 IM1CT=1,NROW
      IM1=IM1CT+LMN
      SGN=-SGN
      CSUM=D(IM1CT,1)
      MULT=1
520   IF(DABS(CSUM).LT.DLBIG) GOTO 510
      MULT=MULT*2
      CSUM=0.5*CSUM
      GOTO 520
510   FAC=10.D0**CSUM
      SFAC=SMALL/FAC
      NM2=MIN0(IM1-1,NM-IM1)
      NM2P1=NM2+1
      DO 130 IM2=1,NM2P1
      IF((D(IM1CT,IM2+1).EQ.0.D0).OR.(IM2.GE.NM2)) GOTO 250
      CSUM=CSUM*DFLOAT(MULT)+D(IM1CT,IM2+1)
      MULT=1
220   IF(DABS(CSUM).LT.DLBIG) GOTO 210
      MULT=MULT*2
      CSUM=0.5D0*CSUM
      GOTO 220
210   FAC=10.D0**CSUM
      SFAC=SMALL/FAC
250   IN2=NMP1-IM2
      DO 270 I=1,MULT
      TEMP=D(IM1CT,IN2)
      RMOD=DABS(TEMP)
      IF(RMOD.GT.SFAC) GOTO 260
      D(IM1CT,IN2)=0.D0
      GOTO 130
260   D(IM1CT,IN2)=D(IM1CT,IN2)*FAC
270   CONTINUE
      D(IM1CT,IN2)=SGN*D(IM1CT,IN2)
130   CONTINUE
120   CONTINUE
C
C       FILL REST OF MATRIX
C
400   IF(ISUP.GT.0) GOTO 410
      SGN=-1.D0
      IF((LMN/2)*2.NE.LMN) SGN=1.D0
      DO 420 IM1CT=1,NROW
      SGN=-SGN
      IM1=IM1CT+LMN
      NM2=MIN0(IM1,NMP1-IM1)
      DO 420 IN2=1,NM2
      IM2=NMP1-IN2
420   D(IM1CT,IN2)=SGN*D(IM1CT,IM2)
      DO 430 IM1CT=1,NROW
      IM1=IM1CT+LMN
      IN1=NMP1-IM1
      IN1CT=IN1-LMN
      SGN=-1.D0
      NM2=MIN0(IM1,IN1)
      DO 440 NIT=1,NM2
      SGN=-SGN
      IM2=1+NM2-NIT
      IN2=NMP1-IM2
      IM2CT=IM2-LMN
      IN2CT=IN2-LMN
      D(IN1CT,IN2)=SGN*D(IM1CT,IM2)
      IF(IN2CT.GT.NROW) GOTO 440
      D(IM2CT,IM1)=D(IN1CT,IN2)
      D(IN2CT,IN1)=D(IM1CT,IM2)
440   CONTINUE
430   CONTINUE
      RETURN
410   DO 450 IM1CT=1,NROW
      IM1=IM1CT+LMN
      IN1=NMP1-IM1
      IN1CT=IN1-LMN
      SGN=-1.D0
      NM2=MIN0(IM1,IN1)
      DO 460 NIT=1,NM2
      SGN=-SGN
      IM2=NM-NM2+NIT
      IN2=NMP1-IM2
      IM2CT=IM2-LMN
      IN2CT=IN2-LMN
      D(IN1CT,IN2)=SGN*D(IM1CT,IM2)
      IF(IM2CT.GT.NROW) GOTO 460
      D(IM2CT,IM1)=D(IN1CT,IN2)
      D(IN2CT,IN1)=D(IM1CT,IM2)
460   CONTINUE
450   CONTINUE
      RETURN
      END
