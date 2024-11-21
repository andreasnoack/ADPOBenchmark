using Pumas, CSV, DataFrames

# Define model
mdl = @model begin
  @param begin
    tvvmax ∈ RealDomain(lower = 100, init = 4000)
    tvkm ∈ RealDomain(lower = 100, init = 4000)
    tvv2 ∈ RealDomain(lower = 0, init = 50)
    tvka ∈ RealDomain(lower = 0, init = 1.2)
    θ5 ∈ RealDomain(lower = 0, init = 1.6)
    θ6 ∈ RealDomain(lower = 0, init = 0.7)
    θ7 ∈ RealDomain(lower = 0, init = 1)
    θ8 ∈ RealDomain(lower = 0.0001, init = 0.1)
    θ9 ∈ RealDomain(lower = 0.0001, init = 0.05)
    Ω ∈ PDiagDomain(init = [0.1, 0.1, 0.1])
    σ_prop ∈ RealDomain(lower = 0, init = 0.3)
    σ_add ∈ RealDomain(lower = 0, init = 1)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @pre begin
    Vmax = tvvmax * exp(η[1])
    Km = tvkm * exp(η[2])
    V2 = tvv2 * exp(η[3])
    Ka = tvka
    K23 = θ6
    K32 = θ7
    K24 = θ8
    K42 = θ9
    SC = V2
  end

  @covariates WT

  @dynamics begin
    Depot' = -Ka * Depot
    Central' = Ka * Depot - (Vmax * (Central / SC)^θ5) / (Km^θ5 + (Central / SC)^θ5) - K23 * Central + K32 * Peri1 - K24 * Central + K42 * Peri2
    Peri1' = K23 * Central - K32 * Peri1
    Peri2' = K24 * Central - K42 * Peri2
  end

  @derived begin
    conc = @. Central / SC
    IOBS ~ @. Normal(conc, sqrt((conc^2 * σ_prop) + σ_add))
  end
end

data_df = CSV.read(
  "data/sim_2_5_0_1.csv",
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
  covariates = [:WT],
  amt = :AMT,
  evid = :EVID2,
  cmt = :CMT
)

# $TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Run1Preds.dat

# ;; Phenotype: ([('COMP', 0), ('ETAs', 0), ('V~WT', 0), ('GAMMA', 0)])
# ;; Genotype: [0, 0, 0, 0]
# ;; Num non-influential tokens: 0

for sub in data_pd
  @info sub.id
  loglikelihood(mdl, sub, init_params(mdl), FOCE())
end

loglikelihood(mdl, data_pd, init_params(mdl), FOCE())

ft_foce = @time fit(mdl, data_pd, init_params(mdl), FOCE(), ensemblealg = EnsembleSerial())
ft_gfoce = @time fit(mdl, data_pd, init_params(mdl), Pumas.GFOCE(), ensemblealg = EnsembleSerial())

tmp = [loglikelihood(mdl, sub, init_params(mdl), Pumas.GFOCE()) for sub in data_pd]

_infer = @time infer(ft)
