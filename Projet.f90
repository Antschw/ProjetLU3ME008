PROGRAM Projet 

  USE Projet_mod 

	IMPLICIT NONE
	
	INTEGER, PARAMETER :: nx=30, tfin=120 ! 120s=2min
	REAL :: dt, dx, M, t, err, eps=1.E-5
	REAL, PARAMETER :: T0=750., Tbain=25., L=0.03, alpha=1.E-5
	REAL, DIMENSION(:), ALLOCATABLE :: b, Ta, Tp,c1, x, y, b1, Tn
	REAL, DIMENSION(:,:), ALLOCATABLE :: A, C, Lu, U
  INTEGER, DIMENSION(8) :: T1,T2,T3
  INTEGER :: i, j, n, nmax, cpt, dmin, ds, dms, TF
  
	
	dt=30./(nx**2)        !Calcul de dt de base
	dx=L/nx               !Calcul de dx
	M=(dx**2)/(alpha*dt)	!Calcul de M
  nmax=NINT(tfin/dt)		!Calcul de nmax
  cpt=0				          !Initialisation compteur
  

	ALLOCATE(b(nx-1), Ta(nx-1), A(nx-1,nx-1), Tp(nx-1))								                                    !Allocation variables Schéma Explicite
  ALLOCATE(C(nx-1,nx-1), b1(nx-1), c1(nx-1), Lu(nx-1,nx-1), U(nx-1,nx-1), x(nx-1), y(nx-1), Tn(nx-1))		!Allocation variables Schéma Implicite

  OPEN (10, FILE='graph_30.dat')	!Ouverture des fichiers
  OPEN (11, FILE='graph_60.dat')
  OPEN (12, FILE='graph_90.dat')
  OPEN (13, FILE='graph_120.dat')
  OPEN (14, FILE='graph_undemi.dat')
  OPEN (15, FILE='graph_profil.dat')
	OPEN (16, FILE='graph_normta.dat')
	OPEN (17, FILE='graph_normx.dat')

  CALL DATE_AND_TIME(VALUES=T1) !Marqueur de temps T1
  
	! Schéma Explicite

  DO 

  	b=0.							                      	!Initialisation à 0 des valeurs de b
  	b(1)=1/M*Tbain ; b(nx-1)=1/M*Tbain 				!Oubli du 1/M dans les conditions initiales
  	A=0.								                      !Initalisation des valeurs de la matrice A

  	DO i=1 ,nx-2     	!Calcul des valeurs de A         
    	A(i,i)=1-(2/M)
  		A(i+1,i)=1/M
    	A(i,i+1)=1/M
    END DO  

  	A(nx-1,nx-1)=1-2/M


    Ta=T0 		!Initialisation T(xi,tn)=T0
    cpt=cpt+1 	!Compteur pour tracer les graphs que pour dt d'origine

    DO n = 1, nmax
    
      Ta=MATMUL(A,ta)+b     !Temperature actuelle

      IF (cpt==1) THEN

        WRITE(14,FMT=*) n*dt, Ta(15) 	!Température x=L/2 en fonction du temps
    
        !!! Question 3 !!!
    
        IF((n*dt)==30) THEN 	!Température dans la plaque à n*dt=30s
          DO i=1,nx-1
           WRITE(10, FMT=*) i*dx,Ta(i)
          END DO
        END IF

       IF((n*dt)==60) THEN 	!Température dans la plaque à n*dt=60s
          DO i=1,nx-1
           WRITE(11, FMT=*) i*dx,Ta(i)
          END DO
        END IF

        IF((n*dt)==90) THEN 	!Température dans la plaque à n*dt=90s
          DO i=1,nx-1
            WRITE(12, FMT=*) i*dx,Ta(i)
          END DO
        END IF

        IF((n*dt)==120) THEN 	!Température dans la plaque à n*dt=120s
       	  DO i=1,nx-1
            WRITE(13, FMT=*) i*dx,Ta(i)
          END DO
        END IF

      END IF  

    END DO 

    err=ABS(SUM(ABS(Ta))-SUM(ABS(Tp)))

	  WRITE (16, FMT=*) dt,SUM(ABS(Ta))
    
    IF (err<eps) EXIT !Condition de sortie de la boucle

    Tp=Ta 					        !Tp sauvegarde anciennes valeurs Ta
	  dt=dt/5 				        !Nouveau dt Q4
    M=(dx**2)/(alpha*dt)		!Nouveau M associé à dt
    nmax=NINT(tfin/dt)			!Nouveau nmax associé à dt

  END DO

  PRINT *, "L'ERREUR POUR EPSILON =", eps, "EST DE : ",err !Affichage de l'erreur pour un epsilon fixé

  CALL DATE_AND_TIME(VALUES=T2) !Marqueur de temps T2

  CALL dift(T1,T2,6,TF) ; dmin=TF !Calcul temps passé entre T1 et T2
  CALL dift(T1,T2,7,TF) ; ds=TF
  CALL diftms(T1,T2,8,TF) ; dms=TF
  PRINT*, "TEMPS EXECUTION SCHEMA EXPLICIT : ",dmin,"min",ds, "s",dms,"ms"

  ! Schéma Implicite

 	DO

    x=0.	!Initialisation des valeurs de x
    C=0.	!Initialisation de la matrice C

    DO i=1 ,nx-2                !Calcul de la matrice C
      C(i,i)=1+(2/M)
      C(i+1,i)=-(1/M)
      C(i,i+1)=-(1/M)
    END DO  

    C(nx-1,nx-1)=1+(2/M)	

    c1=0.					! Calcul des valeurs de c1
    c1(1)=(1/M)*Tbain		
    c1(nx-1)=(1/M)*Tbain	

    Tn=750.					!Initialisation des valeurs de Tn

    b1=Tn+c1				!Initialisation des valeurs de b1

    !Factorisation LU

    DO n = 1, nmax

      x=0.	!Réinitialisaton des valeurs de x
      y=0.	!Initialisaton des valeurs de y
      
      !Décomposition A=LU
      CALL decomp(C,Lu,U)

      !Résolution y=Lc
      CALL descente(Lu,b1,y)

      !Résolution x=Uy
      CALL remontee(U,y,x)
			
			IF (cpt==1) THEN
				WRITE(15,FMT=*) n*dt, x(15) !Température x=L/2 en fonction du temps
			END IF

      Tn=x		  !Sauvegarde de la valeur de Tn
      b1=Tn+c1	!Calcul du nouveau b1
     
    END DO
	
		WRITE(17, FMT=*) dt,SUM(ABS(x))
	

    IF (cpt==1) EXIT		!Condition de sortie de la boucle 

    dt=dt*5 				        !Nouveau dt Q8
    M=(dx**2)/(alpha*dt)		!Changement de M avec dt
    nmax=NINT(tfin/dt)			!Changement de nmax avec dt
    cpt=cpt-1				        !Compteur des valeurs de dt

  END DO
	 
  PRINT*,"RE-OBTENTION DU DT D'ORIGINE :", dt

  CALL DATE_AND_TIME(VALUES=T3) !Marqueur de temps T3

  CALL dift(T2,T3,6,TF) ; dmin=TF !Calcul temps passé entre T2 et T3
  CALL dift(T2,T3,7,TF) ; ds=TF
  CALL diftms(T2,T3,8,TF) ; dms=TF
  PRINT*, "TEMPS EXECUTION SCHEMA IMPLICIT : ",dmin,"min",ds, "s",dms,"ms"

  CLOSE (10); CLOSE (11); CLOSE (12); CLOSE (13); CLOSE(14); CLOSE(15); CLOSE(16); CLOSE(17)	!Fermeture des fichiers

	DEALLOCATE(b,Ta,A,Tp,C,c1,Lu,U,x,y,b1,Tn)	!Désallocation des variables
  
END PROGRAM Projet
