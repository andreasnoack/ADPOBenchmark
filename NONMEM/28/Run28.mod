$SIZES MAXFCN = 5000000
$PROBLEM    Dual Numbers benchmark
$INPUT       ID TIME AMT DV WT
  ;; $DATA      c:/git/adpoBenchmark/data/sim_1_3_1_0.csv IGNORE=@
$DATA      C:/git/adpoBenchmark/data/sim_1_3_1_0.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=7
$MODEL
  COMP=(DEPOT,DEFDOSE)
  COMP=(CENTRAL,NODOSE,DEFOBS)
  COMP=(PERI,NODOSE)
$PK
  CWT = WT/70

  TVVMAX= THETA(1)
  VMAX=TVVMAX*EXP(ETA(1)) 
  TVKM = THETA(2)
  KM = TVKM *EXP(ETA(2)) 
  TVV2=THETA(3) *CWT**THETA(5)
  V2=TVV2 *EXP(ETA(3)) 
  TVKA=THETA(4)
  KA=TVKA *EXP(ETA(4)) 
  K23=THETA(6)
  K32=THETA(7)
  SC = V2
$ERROR
  IPRED = F
  Y=F*EXP(EPS(2)) + EPS(1)
$DES
  CONC = A(2)/SC
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) -K23*A(2)+K32*A(3)
  DADT(3) = K23*A(2)-K32*A(3)

$THETA
  (100,4000)	; THETA(1) VMAX UNITS =  MG/HR
  (100,4000)	; THETA(2) KM UNITS =  CONC
  (0,50) 	; THETA(3) V  UNITS = L
  (0,1.2) 	; THETA(4) KA UNITS = 1/HR

  (0,1.1) 	;; THETA(5) V~WT
  (0.0001,2)	 ;; THETA(6) K23
  (0.0001,3)	 ;; THETA(7) K32
; empty $OMEGA
$OMEGA 
  0.1		; ETA(1) ETA ON VMAX
  0.1		; ETA(2) ETA ON KM
  0.1		; ETA(3) ETA ON V2
  0.1		; ETA(4) ETA ON KA
$SIGMA
  (1) ;; EPS(1) ADDITIVE
  (0.3) ;; EPS(2) PROPORTIONAL
  $EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
  $COV UNCOND PRECOND=2
  ;;; Model Identifier =  1,3,1,0

$TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat

;; Phenotype: ([('COMP', 1), ('ETAs', 3), ('V~WT', 1), ('GAMMA', 0)])
;; Genotype: [1, 3, 1, 0]
;; Num non-influential tokens: 0
