	program GaAs
	
	parameter(NR_STEPS = 10000)
	integer ierr,nivel,nivmax,funconda
	real*8 z,zstep,Sel,z_min
	real*8 massa,potencial,Vx,a0,Ry
	real*8 WFel(NR_STEPS)
	real*8 D(NR_STEPS),E(NR_STEPS)
	real*8 auxiliar
	real*8 ll,ll2,pi,m1,m2,me0,hc2,qel
	REAL*8 ZLPK(NR_STEPS,NR_STEPS)		
	INTEGER IL,IU,IM,ILDZ,INFO
	INTEGER IFAIL(NR_STEPS),IWORK(5*NR_STEPS)
	REAL*8  VL,VU,ABSTOL
	REAL*8  W(NR_STEPS),WORK(5*NR_STEPS)
	CHARACTER JOBZ,RANGE
	
      !NEP  -> nivel de energia do portador    - NEP(nivel)
      !WFel -> função de onda do eletron       - WFel(NR_STEPS)

	parameter (nivmax=5)  !nível de energia máximo que vai ser calculado
      real*8 NEP(nivmax)    
	
	common/massas/m1,m2
	common/limits/ll
	common/AlturaPotenc/Vx
	common/adim/ry,a0
	common/axis/zstep,z_min

      open(unit = 12, file = 'wave.dat')


C       ------- LAPACK ----------
      funconda = 1  ! =0 calcula só energias, =1 calcula também funções de onda
	ABSTOL=0.0d0  ! Tolerância do erro -> deixar 0.0d0
	IL=1       ! IL e IU são as posições dos autoestados
      IU=nivmax
	VL=0.d0    ! VL e VU são os limites do intervalo numérico
 	VU=0.d0
      RANGE='I'  ! =I calcula do IL-ésimo ao IU-ésimo estado
	           ! =V calcula autovalores dentro do intervalo (VL,VU]
C      -------------------------
	   
      Ry = 13.6058D+3 !adimensionaliza as energias  (meV)
      a0 = 0.5292D0   !adimensionaliza as posições (angstron)
      pi = 4.D0 * DATAN( 1.D0 )
	zstep = 0.5d0/a0
	z_min = NR_STEPS*zstep/2. !módulo do valor de z nas extremidades
	m1 = 0.067d0   !massa efetiva na região do poço
	m2 = 0.063d0   !massa efetiva na região da barreira
	ll = 350.d0/a0    !largura do poço quântico
	Vx = 250.d0/Ry    !altura do poço quântico

c       <><><><><><><><><><><><><> PARTE PRINCIPAL <><><><><><><><><><><><>
c       <><><><><>< CALCULO DA FUNCAO DE ONDA E NIVEIS DE ENERGIA <><><>><>
	
	do i=1,NR_STEPS,1	!LOOP VARREDURA NO EIXO Z	 
	   
	 z=-z_min + zstep*(i)
	 D(i)=(1./massa(z-zstep) + 2./massa(z) +
     &  	   1./massa(z+zstep))/(2.*zstep*zstep) + potencial(z)
	if(i.ne.1)E(i-1)=-(1./massa(z-zstep)+1./massa(z))/(2.*zstep*zstep)

	enddo		!FIM DO LOOP DA VARREDURA DO EIXO Z
	   
	   
	if(funconda.eq.1)then
	 
	 JOBZ='V'
       write(*,*)"--- calculando a funcao de onda e a energia ---"
c      Input:  Vetores D e E, 
	 call DSTEVX(JOBZ,RANGE,NR_STEPS,D,E,VL,VU,IL,IU,ABSTOL,
     &    IM,W,ZLPK,NR_STEPS,WORK,IWORK,IFAIL,INFO)

c      Output: Vetor W dá os autovalores, ZLPK dá os autovetores
	      
	  nivel=1 !nível = nível de energia q se deseja a função
	          !armazenando as funções de onda
	  do i=1,NR_STEPS
         WFel(i)=ZLPK(i,nivel)/zstep
	  enddo
		 
	 !encontrando o fator que normaliza
	  call Norma(NR_STEPS,WFel,Sel)
		      
	      
	 do i=1,NR_STEPS !dividindo pelo fator, para normalizar
	   z=-z_min + zstep*(i+1d-8)
 	   write(12,20)z,WFel(i),(WFel(i)/(Sel))**2
	   WFel(i)=WFel(i)/(Sel)
	 enddo
	      
 20	 format(3e16.6)
	      
	else  !else de "if funconda eq 1"
	     
       JOBZ='N'
       write(*,*)"--- calculando somente e energia ---"
c      Input:  Vetores D e E                 
	 call DSTEVX(JOBZ,RANGE,NR_STEPS,D,E,VL,VU,IL,IU,ABSTOL,
     &         IM,W,ZLPK,NR_STEPS,WORK,IWORK,IFAIL,INFO)
c      Output: Vetor W dá os autovalores
	      
	endif

	auxiliar=potencial(z_min)*Ry  !deixa o potencial no ultimo
	                                    !ponto da barreira em meV
	   
	do nivel=1,nivmax
	 NEP(nivel)=W(nivel)*Ry !deixa a energia em meV  
	 if(NEP(nivel).gt.auxiliar)then !comparando o nível com o potencial da barreira
	  NEP(nivel)=0.0 !se o nível for mais alto q a barreira, ele não conta, por isso é posto como zero
	 endif
      enddo
	   
c       <><><><><>< FIM DA PARTE PRINCIPAL <><><><><><><><><><><><><><>
  
c     ------------ ESCREVENDO OS AUTOESTADOS NA TELA ------------------	   
	do nivel=1,nivmax	
	 write(*,*)nivel,NEP(nivel) 
        enddo			
!        write(*,*)'o que deveriam ser:',Ry*1.D0*pi**2/ll**2, 
!     &              Ry*4.D0*pi**2/ll**2, Ry*9.D0*pi**2/ll**2 
	
 100	end program
	

C     =================================================================
C                FUNÇÕES QUE DESCREVEM POTENCIAL E MASSA 
C     =================================================================	
	real*8 function massa(z)
	real*8 z
	real*8 m1,m2,ll
	common/massas/m1,m2
	common/limits/ll
	
	if(dabs(z).le.ll/2.)then
	 massa=m1
	else
	 massa=m2
	endif

	end

c     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

	real*8 function potencial(z)
	real*8  z,ry,a0,ll
	real*8  Vx,x
	common/AlturaPotenc/Vx
	common/adim/ry,a0
	common/limits/ll
			
	!---- poço quadrado ----------------
       if(dabs(z).le.ll/2.)then
	  potencial= 0.d0
	 else
	  potencial = Vx 
	 endif 
        
	 !potencial = potencial + q*F*z

    	end

c     <><><><><><><><><><><><><><><><><><><><><><><><><><><>
	subroutine Norma(nx,W,S) !encontra o parâmetro de normalização das funções de onda
	real*8 a,b,h,S,S0,S1
	integer nx,i
	real*8 zstep,z_min,W(nx)    
	
	common/axis/zstep,z_min
	
c	nx=NR_STEPS
	a=-z_min
	b=-z_min+zstep*(nx)
	h=zstep
	S=0.0
	S1=0.0 
	S0=0.0
	
	do i=1,nx
         S0=h*W(i)*W(i)
	 S1=S1+S0
	enddo
	
	S=dsqrt(S1)

	return
	end
