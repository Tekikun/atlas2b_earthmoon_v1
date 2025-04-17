C PROGRAM ATLAS2BGENERAL.F - TABARE GALLARDO 2015
C CALCULATES AN ATLAS OF RESONANCE'S STRENGTHS FOR A GIVEN PLANETARY
C SYSTEM
C www.fisica.edu.uy/~gallardo/atlas/
C gallardo@fisica.edu.uy
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SE(10),SM(10)
      DIMENSION NPLA(100000),NPQ(100000),NP(100000),SEMI(100000)
      DIMENSION A2(100000),NCPLA(100000),NCPQ(100000),NCP(100000)
      CHARACTER*56 cabezal
      CHARACTER*21 cabezal2
      CHARACTER*93 cabezal3
      COMMON SE,SM
      cabezal=' pla  |p+q|:|p|     a(au)      e     i     argper       '
      cabezal2=' <R>      strength'
      cabezal3=cabezal//cabezal2

      OPEN(1,FILE="massivebodies.inp")
      read(1,1004) APESOS
C MASS OF THE STAR IN SOLAR MASSES
      READ(1,*) xm0
      read(1,1004) APESOS
C NUMBER OF PLANETS
      READ(1,*) NTOT
C SEMIMAJOR AXES AND MASSES
      read(1,1004) APESOS
      DO I=1,NTOT
        READ(1,*) SE(I),SM(I)
      ENDDO
      CLOSE(1)

      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)' Given a planetary system this program calculates all '
      WRITE(*,*)' TWO BODY RESONANCES of a massless particle in a '
      WRITE(*,*)' given range of semimajor axes with all the planets '
      WRITE(*,*)' taken one by one. Gallardo (2006, Icarus 184, 29).'
      WRITE(*,*)'          www.fisica.edu.uy/~gallardo/atlas'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'  '
      WRITE(*,*)'Put the data for the central star and planets in the '
      WRITE(*,*)'input file              massivebodies.inp  '
      WRITE(*,*)'  '
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)' '

C MAX ORDER
      WRITE(*,*)'max order q?'
      READ(*,*)IQF
C MAX DEGREE
      WRITE(*,*)'max degree |p|?'
      READ(*,*)IPF
C MAX |P+Q|
      WRITE(*,*)'max |p+q|?'
      READ(*,*)NLIMPQ


      
C PLANETS
      WRITE(*,*)'initial and final planet? (ex: 2,7) '
      READ(*,*) IPLA1,IPLA2

      
C LIMITS IN SEMIMAJOR AXIS
      WRITE(*,*)'interval in au: a_min,a_max? (ex: 3.1,3.3)'
      READ(*,*)A2MIN,A2MAX
C ORBITAL ELEMENTS OF THE PARTICLE
      WRITE(*,*)'e,i,argper for the particle? (ex: 0.3,35,127)'
      READ(*,*)exc,ync,arg
C PRECISION IN CALCULATION OF SR
      isimax=36
C MASS OF THE PARTICLE =0
      XM2=0.0D0
      ZM2=XM0+XM2
C COUNTER I
      I=0
C FOR EACH PLANET.....
      DO 200 JPLA=IPLA1,IPLA2
      XM1=SM(JPLA)
      ZM1=XM0+XM1

      A1=SE(JPLA)
C FOR EACH ORDER...
      DO IQ=0,IQF
        Q=IQ*1.D0
C FOR EACH DEGREE FOR INTERIOR RESONANCES....
        DO IP=1,IPF
          P=IP*1.D0
C LIMIT |P+Q|
          isum=abs(iq+ip)
          if(isum.le.nlimpq) then
C CALCULATE SEMIMAJOR AXIS
          sem=((P/(P+Q))**2*ZM2/ZM1)**(1.D0/3.D0)*A1
            IF (sem.LE.A2MAX.AND.sem.GE.A2MIN)  THEN
              i=i+1
              NPLA(I)=JPLA
              NPQ(I)=INT(ABS(P+Q))
              NP(I)=INT(ABS(P))
              A2(I)=sem
            endif
          endif
        ENDDO
      ENDDO
C NOW THE SAME FOR EXTERIOR RESONANCES
      DO IQ=0,IQF
        Q=IQ*1.D0
        DO IP=IQ,IPF-1
          P=-(IP+1)*1.D0
C LIMIT |P+Q|
          isum=abs(iq-ip-1)
          if(isum.le.nlimpq) then
          sem=((P/(P+Q))**2*ZM2/ZM1)**(1.D0/3.D0)*A1
            IF (sem.LE.A2MAX.AND.sem.GE.A2MIN)  THEN
              I=I+1
              NPLA(I)=JPLA
              NPQ(I)=INT(ABS(P+Q))
              NP(I)=INT(ABS(P))
              A2(I)=sem
            endif
          endif
        ENDDO
      ENDDO
 200  CONTINUE
C NUMBER OF POSSIBLE (redundant) RESONANCES, MUST BE LESS THAN 100000
      IDAT=I

      if(idat.ge.100000) then
        write(*,*)'error: too much resonances'
        stop
      end if
C ELIMINATE REPEATED, EXAMPLE: 4:2 IS THE SAME AS 2:1
      DO I=1,IDAT
        DO J=I+1,IDAT
          IF(NPLA(J).EQ.NPLA(I)) THEN
            IF (A2(J).EQ.A2(I).AND.NPQ(J).GE.NPQ(I)) THEN
              A2(J)=0.D0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      NCONT=0
      DO I=1,IDAT
        IF (A2(I).NE.0.D0) THEN
C LIMIT |P+Q|
C          IF (NPQ(I).LE.NLIMPQ)  THEN
            NCONT=NCONT+1
            NCPQ(NCONT)=NPQ(I)
            NCP(NCONT)=NP(I)
            NCPLA(NCONT)=NPLA(I)
            SEMI(NCONT)=A2(I)
C          ENDIF
        ENDIF
      ENDDO
C RESONANCES REPEATED WERE ELIMINATED
      WRITE(*,*)'I WILL CALCULATE ',NCONT,' RESONANCES'
C OK. WE NEED TO CALCULATE THE STRENGTH OF NCONT RESONANCES
      OPEN(1,FILE='atlas2bres.dat',STATUS='UNKNOWN',ACCESS='APPEND')
      write(1,*)"strength = semiamplitude of R in solar mass,au,day"
      write(1,*) cabezal3
      DO I=1,NCONT
        VPQ=DFLOAT(NCPQ(I))
        VP=DFLOAT(NCP(I))
        CALL FUERZA(ISIMAX,VPQ,VP,NCPLA(I),EXC,YNC,ARG,FME,FUE)
        WRITE(1,45) NCPLA(I),NCPQ(I),NCP(I),SEMI(I),EXC,YNC,ARG,FME,FUE
        WRITE(*,45) NCPLA(I),NCPQ(I),NCP(I),SEMI(I),EXC,YNC,ARG,FME,FUE
      ENDDO
      CLOSE(1)
 45   FORMAT(3I5,F13.6,3F7.2,2(1PE14.4))
1004  FORMAT(A76)
      STOP
      END
C --------------------------------------------------------------------
C CALCULATES THE STRENGTH OF R(sigma).
C PLANETS ASSUMED IN CIRCULAR AND ZERO INCLINATION ORBITS
C THIS SUBROUTINE IS A SLIGTH MODIFICATION OF RSIGMA.F
      SUBROUTINE FUERZA(ISIMAX,NPQ,NP,INPLA,EXC,YNC,ARGPER,FME,FUE)
      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION RP(400),SE(10),SM(10)
      COMMON SE,SM
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      ERROR=1.D-12
      KGAUS=0.01720209895D0
      KG2=KGAUS**2
C NPQ=|P+Q|
C NP=|P|
C ORDER q
      ORQ=DABS(NPQ-NP)
      MAXFAC=NPQ
      IF(NPQ.LE.NP) THEN
C IT IS AN EXTERIOR RESONANCE  OR TROJANS
        MAXFAC=NP
        NPQ=-NPQ
        NP=-NP
      ENDIF
C sigma is invariant to changes in aries then I take aries at the node
C so longnode(particle)=0
C The planet is in circular and i=0 orbit  then the longper is not defined
C so I take \omega=\Omega=0 for the planet
      LONOD=0.D0
      LOPER=argper
C THE PLANET IS 1
      A1=SE(INPLA)
      E1=0.D0
      J1G=0.D0
      L1G=0.D0
      P1G=0.D0
      J1=J1G*G2R
      L1=L1G*G2R
      P1=P1G*G2R
C PARTICLE IS 2
      A2=A1/(DABS(NPQ/NP))**(2.D0/3.D0)
      E2=EXC
      J2GRA=YNC
      J2=J2GRA*G2R
      L2G=LONOD
      P2G=LOPER
      L2=L2G*G2R
      P2=P2G*G2R
C NUMBER OF EVALUATIONS OF R(sigma) BETWEEN O AND 360
C FOR MORE PRECISION USE MORE POINTS:
C 18 IS ENOUGHT FOR A ROUGH ESTIMATION
C 180 IS EXAGGERATEDLY PRECISE
C      ISIMAX=18
C      ISIMAX=180
C STEPS OF THE NUMERICAL INTEGRATION
C FOR MORE PRECISTION USE MORE POINTS
C TAKE CARE!!!!:     THIS INTEGRATION *MUST* BE PRECISE
      IPASOS=100*INT(MAXFAC)
C sigma = (p+q)*lambda_p - p*lambda - q*varpi
      DO ISI=1,ISIMAX
C acrit is sigma
        ACRIT=DFLOAT(ISI-1)/DFLOAT(ISIMAX)*360.D0
        TETA=ACRIT+ORQ*LOPER
 901    IF (TETA.GT.360.D0) THEN
          TETA=TETA-360.D0
          GOTO 901
        END IF
C to radians
        TETAR=TETA*G2R
        RTOT=CERO
        DO ILAMBDA1=1,IPASOS
C all possible configurations occur after p revolutions of the planet
C LAMBDA PLANET
        LA1=DFLOAT(ILAMBDA1)/DFLOAT(IPASOS)*TWOPI*DABS(NP)
C MEAN ANOMALY
        AM1=LA1-P1
 200    IF (AM1.GT.TWOPI) THEN
          AM1=AM1-TWOPI
          GOTO 200
        END IF
 201    IF (AM1.LT.CERO) THEN
          AM1=AM1+TWOPI
          GOTO 201
        END IF
C GIVEN LAMBDA PLANET CALCULATE LAMBDA PARTICLE
        LA2=(NPQ*LA1-TETAR)/NP
C MEAN ANOMALY PARTICLE
        AM2=LA2-P2
 202    IF (AM2.GT.TWOPI) THEN
          AM2=AM2-TWOPI
          GOTO 202
        END IF
 203    IF (AM2.LT.CERO) THEN
          AM2=AM2+TWOPI
          GOTO 203
        END IF
C CRUDE METHOD FOR SOLVING KEPLER FOR PLANET (A WASTE OF TIME.....)
        AEX=AM1
 210    AEXN=AM1+E1*DSIN(AEX)
        IF(DABS(AEXN-AEX).GT.ERROR) THEN
          AEX=AEXN
          GOTO 210
        ENDIF
        AEX1=AEXN
C TRUE ANOMALY
        AVE1=DACOS((DCOS(AEX1)-E1)/(UNO-E1*DCOS(AEX1)))
        IF(AM1.GT.PI) AVE1=TWOPI-AVE1
C HELIOCENTRIC DISTANCE
        R1=A1*(UNO-E1*DCOS(AEX1))
C CRUDE METHOD FOR SOLVING KEPLER FOR PARTICLE (ECCENTRIC)
        AEX=AM2
 211    AEXN=AM2+E2*DSIN(AEX)
        IF(DABS(AEXN-AEX).GT.ERROR) THEN
          AEX=AEXN
          GOTO 211
        ENDIF
        AEX2=AEXN
C TRUE ANOMALY
        AVE2=DACOS((DCOS(AEX2)-E2)/(UNO-E2*DCOS(AEX2)))
        IF(AM2.GT.PI) AVE2=TWOPI-AVE2
C HELIOCENTRIC DISTANCE
        R2=A2*(UNO-E2*DCOS(AEX2))
C RADIUS VECTOR FOR PLANET
        X1=R1*(DCOS(L1)*DCOS(P1-L1+AVE1)-DSIN(L1)*DSIN(P1-L1+AVE1)
     **DCOS(J1))
        Y1=R1*(DSIN(L1)*DCOS(P1-L1+AVE1)+DCOS(L1)*DSIN(P1-L1+AVE1)
     **DCOS(J1))
        Z1=R1*DSIN(P1-L1+AVE1)*DSIN(J1)
C RADIUS VECTOR FOR PARTICLE
        X2=R2*(DCOS(L2)*DCOS(P2-L2+AVE2)-DSIN(L2)*DSIN(P2-L2+AVE2)
     **DCOS(J2))
        Y2=R2*(DSIN(L2)*DCOS(P2-L2+AVE2)+DCOS(L2)*DSIN(P2-L2+AVE2)
     **DCOS(J2))
        Z2=R2*DSIN(P2-L2+AVE2)*DSIN(J2)
C DISTANCE
        DELTA=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
C SCALAR PRODUCT R1*R2
        R1R2=X1*X2+Y1*Y2+Z1*Z2
C DIRECT PART
        RDIR=UNO/DELTA

C      if(rdir.gt.100.d0) rdir=0.d0



C INDIRECT PART (PROBABLY THIS IS A WASTE OF TIME)
        RIND=-R1R2/R1**3
C DISTURBING FUNCTION
        RPER=RDIR+RIND
C SUMATION FOR A CRUDE INTEGRAL
        RTOT=RTOT+RPER
        ENDDO
C END OF LOOP FOR LAMBDA
      RTOT=RTOT/DFLOAT(IPASOS)
C NOW MULTIPLY BY MASS AND CONSTANT OF GRAVITATION
      RP(ISI)=RTOT*SM(INPLA)*kg2
      ENDDO
C END OF LOOP FOR SIGMA
C CALCULATION OF MEAN VALUE <R>
      VMEDIO=0.D0
      DO I2=1,ISIMAX
        VMEDIO=VMEDIO+RP(I2)
      ENDDO
      VMEDIO=VMEDIO/DFLOAT(ISIMAX)
C FIND R_max AND R_min
      VRMAX=-999999.D0
      VRMIN=999999.D0
      DO I2=1,ISIMAX
        IF(RP(I2).LT.VRMIN) THEN
          VRMIN=RP(I2)
        ENDIF
        IF(RP(I2).GT.VRMAX) THEN
          VRMAX=RP(I2)
        ENDIF
      ENDDO
      FME=VMEDIO
      FUE=VMEDIO-VRMIN
      RETURN
      END

