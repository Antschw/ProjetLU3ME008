MODULE Projet_mod
  
  IMPLICIT NONE
  
  CONTAINS

  ! Affichage d'une matrice A de dimensions nl x nc
  SUBROUTINE affiche (M)
    
    REAL, DIMENSION (:, :), INTENT(IN)	:: M
    INTEGER 		:: i, j, nl, nc 
    
    nl = SIZE(M,1)
    nc = SIZE(M,2)  
    
    DO i=1, nl
      WRITE (*,fmt='(30E12.4)') (M(i,j), j=1, nc) 
    END DO
    
    PRINT *
  
  END SUBROUTINE affiche

  ! Affichage sur une ligne d'un vecteur x de dimension n
  SUBROUTINE affiche_vec (x)
    
    REAL, DIMENSION (:), INTENT(IN)	:: x
    INTEGER 		:: i,n
    
    n = SIZE(x)
    
    WRITE (*,fmt='(30E12.4)') (x(i), i=1, n)  
    PRINT *
  
  END SUBROUTINE affiche_vec

  !Décomposition d'un systeme en A=LU
  SUBROUTINE decomp(A,L,U)

    REAL, DIMENSION(:,:), INTENT(IN)  :: A
    REAL, DIMENSION(:,:), INTENT(OUT) :: L, U
    REAL :: s
    INTEGER :: i, j, k, n

    n = SIZE(A,1)
    
    DO j=1,n
      L(j,j)=1.
      DO i=1,j
        s=0.
        DO k=1,i-1
          s=s+L(i,k)*U(k,j)
        ENDDO
        U(i,j)=A(i,j)-s
      ENDDO
      DO i=j+1,n
        s=0.
        DO k=1,j-1
          s=s+L(i,k)*U(k,j)
        ENDDO
        L(i,j)=(A(i,j)-s)/U(j,j)
      ENDDO
    ENDDO

  END SUBROUTINE decomp

  ! Resolution d'un systeme triangulaire inférieur
  SUBROUTINE descente(L,b,y)

    REAL, DIMENSION(:,:), INTENT(IN) :: L
    REAL, DIMENSION(:), INTENT(IN)   :: b
    REAL, DIMENSION(:), INTENT(OUT)  :: y
    REAL  	:: s
    INTEGER 	:: i, j, k, n
    
    n = SIZE(L,1)
    
    DO i=1,n 
      s=0.
      DO j=1,i-1
        s=s+L(i,j)*y(j)
      ENDDO
      y(i)=(b(i)-s)/L(i,i)
    ENDDO

  END SUBROUTINE descente

  ! Resolution d'un systeme triangulaire supérieur
  SUBROUTINE remontee(U,y,x)

    REAL, DIMENSION(:,:), INTENT(IN) :: U
    REAL, DIMENSION(:), INTENT(IN)   :: y
    REAL, DIMENSION(:), INTENT(OUT)  :: x
    REAL        :: s
    INTEGER     :: i, j, k, n
    
    n = SIZE(U,1)
    
    DO i=n,1,-1
      s=0.
      DO j=i+1,n
        s=s+U(i,j)*x(j)
      ENDDO
      x(i)=(y(i)-s)/U(i,i)
    ENDDO 

  END SUBROUTINE remontee

  ! Resolution de la difference de temps pour min et s
  SUBROUTINE dift(T1,T2,n,TF)

	  INTEGER, DIMENSION(8), INTENT(IN) :: T1,T2
	  INTEGER, INTENT(IN):: n
	  INTEGER, INTENT(INOUT):: TF

	  IF (T2(n)>=T1(n)) THEN
  	  TF= T2(n) - T1(n)
	  ELSE
  	  TF=(60-T1(n)) + T2(n)
	  END IF

  END SUBROUTINE dift

  ! Resolution de la difference de temps pour ms
  SUBROUTINE diftms(T1,T2,n,TF)

	  INTEGER, DIMENSION(8), INTENT(IN) :: T1,T2
	  INTEGER, INTENT(IN):: n
	  INTEGER, INTENT(INOUT):: TF

	  IF (T2(n)>=T1(n)) THEN
  	  TF= T2(n) - T1(n)
	  ELSE
  	  TF=(1000-T1(n)) + T2(n)
	  END IF

  END SUBROUTINE diftms

END MODULE Projet_mod