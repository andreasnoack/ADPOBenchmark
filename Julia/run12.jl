using Pumas, CSV, DataFrames

# $SIZES MAXFCN = 5000000
# $PROBLEM    Dual Numbers benchmark
# $INPUT       ID TIME AMT DV WT
#   ;; $DATA      c:/git/adpoBenchmark/data/sim_1_5_0_0.csv IGNORE=@
# $DATA      ..\\..\\data\\sim_1_5_0_0.csv IGNORE=@

mdl_react = @model begin

# $SUBROUTINE ADVAN6 TOL=7
# $MODEL
#   COMP=(DEPOT,DEFDOSE)
#   COMP=(CENTRAL,NODOSE,DEFOBS)
#   COMP=(PERI,NODOSE)

# $THETA
#   (100,4000)    ; THETA(1) VMAX UNITS =  MG/HR
#   (100,4000)    ; THETA(2) KM UNITS =  CONC
#   (0,50)    ; THETA(3) V  UNITS = L
#   (0,1.2)   ; THETA(4) KA UNITS = 1/HR


#   (0.0001,2)     ;; THETA(5) K23
#   (0.0001,3)     ;; THETA(6) K32
# ; empty $OMEGA
# $OMEGA BLOCK(3)
#   0.1       ; ETA(1) ETA ON VMAX
#   0.05 0.1      ; ETA(2) ETA ON KM
#   0.05 0.05 0.1     ; ETA(3) ETA ON V2
# $SIGMA
#   (1) ;; EPS(1) ADDITIVE
#   (0.3) ;; EPS(2) PROPORTIONAL
#   $EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
#   $COV UNCOND PRECOND=2
#   ;;; Model Identifier =  1,5,0,0

  @param begin
    θvmax ∈ RealDomain(lower = 100.0, init = 4000.0)
    θkm ∈ RealDomain(lower = 100.0, init = 4000.0)
    θV ∈ RealDomain(lower = 0.0, init = 50.0)
    θka ∈ RealDomain(lower = 0.0, init = 1.2)
    θk23 ∈ RealDomain(lower = 0.0001, init = 2.0)
    θk32 ∈ RealDomain(lower = 0.0001, init = 3.0)

    Ω ∈ PSDDomain(
      init = [
        0.1  0.05 0.05
        0.05 0.1  0.05
        0.05 0.05 0.1
      ]
    )

    σₐ ∈ RealDomain(lower = 0.0, init = 1.0)
    σₚ ∈ RealDomain(lower = 0.0, init = sqrt(0.3))
  end

# $PK
#   CWT = WT/70

#   TVVMAX= THETA(1)
#   VMAX=TVVMAX*EXP(ETA(1))
#   TVKM = THETA(2)
#   KM = TVKM *EXP(ETA(2))
#   TVV2=THETA(3)
#   V2=TVV2 *EXP(ETA(3))
#   TVKA=THETA(4)
#   KA=TVKA
#   K23=THETA(5)
#   K32=THETA(6)
#   SC = V2

  @random begin
    η ~ MvLogNormal([
      log(θvmax),
      log(θkm),
      log(θV)
      ],
      Ω
    )
  end

  @pre begin
    vmax = η[1]
    km = η[2]
    V = η[3]
  end

# $DES
#   CONC = A(2)/SC
#   DADT(1) = -KA*A(1)
#   DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) -K23*A(2)+K32*A(3)
#   DADT(3) = K23*A(2)-K32*A(3)

  @vars begin
    conc = Central / V
  end

  @reactions begin
    θka, Depot --> Central
    (θk23, θk32), Central <--> Peripheral
    vmax / V / (km + Central/V), Central --> 0
  end

  # @dynamics begin
  #   Depot' = -θka*Depot
  #   Central' = θka*Depot - vmax*conc/(km + conc) - k23*Central + k32*Peripheral
  #   Peripheral' = k23*Central - k32*Central
  # end

# $ERROR
#   IPRED = F
#   Y=F*EXP(EPS(2)) + EPS(1)
  @derived begin
    IOBS ~ @. Normal(conc, sqrt(σₐ^2 + (conc*σₚ)^2))
  end
end

mdl_ode = @model begin

# $SUBROUTINE ADVAN6 TOL=7
# $MODEL
#   COMP=(DEPOT,DEFDOSE)
#   COMP=(CENTRAL,NODOSE,DEFOBS)
#   COMP=(PERI,NODOSE)

# $THETA
#   (100,4000)    ; THETA(1) VMAX UNITS =  MG/HR
#   (100,4000)    ; THETA(2) KM UNITS =  CONC
#   (0,50)    ; THETA(3) V  UNITS = L
#   (0,1.2)   ; THETA(4) KA UNITS = 1/HR


#   (0.0001,2)     ;; THETA(5) K23
#   (0.0001,3)     ;; THETA(6) K32
# ; empty $OMEGA
# $OMEGA BLOCK(3)
#   0.1       ; ETA(1) ETA ON VMAX
#   0.05 0.1      ; ETA(2) ETA ON KM
#   0.05 0.05 0.1     ; ETA(3) ETA ON V2
# $SIGMA
#   (1) ;; EPS(1) ADDITIVE
#   (0.3) ;; EPS(2) PROPORTIONAL
#   $EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
#   $COV UNCOND PRECOND=2
#   ;;; Model Identifier =  1,5,0,0

  @param begin
    θvmax ∈ RealDomain(lower = 100.0, init = 4000.0)
    θkm ∈ RealDomain(lower = 100.0, init = 4000.0)
    θV ∈ RealDomain(lower = 0.0, init = 50.0)
    θka ∈ RealDomain(lower = 0.0, init = 1.2)
    θk23 ∈ RealDomain(lower = 0.0001, init = 2.0)
    θk32 ∈ RealDomain(lower = 0.0001, init = 3.0)

    Ω ∈ PSDDomain(
      init = [
        0.1  0.05 0.05
        0.05 0.1  0.05
        0.05 0.05 0.1
      ]
    )

    σₐ ∈ RealDomain(lower = 0.0, init = 1.0)
    σₚ ∈ RealDomain(lower = 0.0, init = sqrt(0.3))
  end

# $PK
#   CWT = WT/70

#   TVVMAX= THETA(1)
#   VMAX=TVVMAX*EXP(ETA(1))
#   TVKM = THETA(2)
#   KM = TVKM *EXP(ETA(2))
#   TVV2=THETA(3)
#   V2=TVV2 *EXP(ETA(3))
#   TVKA=THETA(4)
#   KA=TVKA
#   K23=THETA(5)
#   K32=THETA(6)
#   SC = V2

  @random begin
    η ~ MvLogNormal([
      log(θvmax),
      log(θkm),
      log(θV)
      ],
      Ω
    )
  end

  @pre begin
    vmax = η[1]
    km = η[2]
    V = η[3]
  end

# $DES
#   CONC = A(2)/SC
#   DADT(1) = -KA*A(1)
#   DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC) -K23*A(2)+K32*A(3)
#   DADT(3) = K23*A(2)-K32*A(3)

  @vars begin
    conc = Central / V
  end

  # @reactions begin
  #   θka, Depot --> Central
  #   (θk23, θk32), Central <--> Peripheral
  #   vmax / V / (km + Central/V), Central --> 0
  # end

  @dynamics begin
    Depot' = -θka*Depot
    Central' = θka*Depot - vmax*conc/(km + conc) - θk23*Central + θk32*Peripheral
    Peripheral' = θk23*Central - θk32*Peripheral
  end

# $ERROR
#   IPRED = F
#   Y=F*EXP(EPS(2)) + EPS(1)
  @derived begin
    IOBS ~ @. Normal(conc, sqrt(σₐ^2 + (conc*σₚ)^2))
  end
end


# $TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat

# ;; Phenotype: ([('COMP', 1), ('ETAs', 5), ('V~WT', 0), ('GAMMA', 0)])
# ;; Genotype: [1, 5, 0, 0]
# ;; Num non-influential tokens: 0

data_df = CSV.read(
  "data/sim_1_5_0_0.csv",
  DataFrame,
  missingstring = ["."]
)
transform!(
  data_df,
  "ID" => ByRow(t -> 1) => "CMT",
  "EVID" => ByRow(t -> t == 2 ? 0 : 1) => "EVID2"
)

data_pd = read_pumas(
  data_df,
  id = :ID,
  time = :TIME,
  observations = [:IOBS],
  amt = :AMT,
  evid = :EVID2,
  cmt = :CMT
)

loglikelihood(mdl_react, data_pd, init_params(mdl), FOCE())
loglikelihood(mdl_ode, data_pd, init_params(mdl), FOCE())

ft_foce = @time fit(mdl_ode, data_pd, init_params(mdl), FOCE(), ensemblealg = EnsembleSerial())
ft_gfoce = @time fit(mdl, data_pd, init_params(mdl), Pumas.GFOCE(), ensemblealg = EnsembleSerial())

[loglikelihood(mdl, sub, init_params(mdl), Pumas.GFOCE()) for sub in data_pd]

_infer = @time infer(ft)
