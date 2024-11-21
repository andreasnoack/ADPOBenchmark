using Pumas, CSV, DataFrames

# $SIZES MAXFCN = 5000000
# $PROBLEM    Dual Numbers benchmark
# $INPUT       ID TIME AMT DV WT
  # ;; $DATA      c:/git/adpoBenchmark/data/sim_0_0_0_0.csv IGNORE=@
# $DATA     ..\\..\\data\\sim_0_0_0_0.csv IGNORE=@

# $SUBROUTINE ADVAN6 TOL=7
# $MODEL
  # COMP=(DEPOT,DEFDOSE)
  # COMP=(CENTRAL,NODOSE,DEFOBS)

mdl = @model begin

# $THETA
#   (100,4000)    ; THETA(1) VMAX UNITS =  MG/HR
#   (100,4000)    ; THETA(2) KM UNITS =  CONC
#   (0,50)    ; THETA(3) V  UNITS = L
#   (0,1.2)   ; THETA(4) KA UNITS = 1/HR
# ; empty $OMEGA
# $OMEGA
#   0.1       ; ETA(1) ETA ON VMAX
# $SIGMA
#   (1) ;; EPS(1) ADDITIVE
#   (0.3) ;; EPS(2) PROPORTIONAL
#   $EST METHOD=COND INTER MAXEVALS=9999 NOHABORT NOOMEGABOUNDTEST NOSIGMABOUNDTEST NOTHETABOUNDTEST
#   $COV UNCOND PRECOND=2
#   ;;; Model Identifier =  0,0,0,0
  @param begin
    θvmax ∈ RealDomain(lower = 100.0, init = 4000.0)
    θkm ∈ RealDomain(lower = 100.0, init = 4000.0)
    θV2 ∈ RealDomain(lower = 0.0, init = 50.0)
    θka ∈ RealDomain(lower = 0.0, init = 1.2)

    ωvmax ∈ RealDomain(lower = 0.0, init = sqrt(0.1))

    σₐ ∈ RealDomain(lower = 0.0, init = 1.0)
    σₚ ∈ RealDomain(lower = 0.0, init = sqrt(0.3))
  end


# $PK
  # CWT = WT/70

  # TVVMAX= THETA(1)
  # VMAX=TVVMAX*EXP(ETA(1))
  # TVKM = THETA(2)
  # KM = TVKM
  # TVV2=THETA(3)
  # V2=TVV2
  # TVKA=THETA(4)
  # KA=TVKA

  # SC = V2
  @random begin
    vmax ~ LogNormal(log(θvmax), ωvmax)
  end

  @pre begin
    km = θkm
    V2 = θV2
    ka = θka
  end

# $DES
  # CONC = A(2)/SC
  # DADT(1) = -KA*A(1)
  # DADT(2) =  KA*A(1)-VMAX*CONC/(KM+CONC)

  @vars begin
    conc = Central / V2
  end

  @dynamics begin
    Depot' = -ka*Depot
    Central' = ka*Depot - vmax * conc / (km + conc)
  end

# $ERROR
  # IPRED = F
  # Y=F*EXP(EPS(2)) + EPS(1)
  @derived begin
    IOBS ~ @. Normal(conc, sqrt(σₐ^2 + (conc*σₚ)^2))
  end
end

data_df = CSV.read(
  "data/sim_0_0_0_0.csv",
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

# $TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat

# ;; Phenotype: ([('COMP', 0), ('ETAs', 0), ('V~WT', 0), ('GAMMA', 0)])
# ;; Genotype: [0, 0, 0, 0]
# ;; Num non-influential tokens: 0

loglikelihood(mdl, data_pd, init_params(mdl), FOCE())
loglikelihood(mdl, data_pd, init_params(mdl), Pumas.GFOCE())

ft_foce = @time fit(mdl, data_pd, init_params(mdl), FOCE(), ensemblealg = EnsembleSerial())
ft_gfoce = @time fit(mdl, data_pd, init_params(mdl), Pumas.GFOCE(), ensemblealg = EnsembleSerial())

tmp = [loglikelihood(mdl, sub, init_params(mdl), Pumas.GFOCE()) for sub in data_pd]

_infer = @time infer(ft)
