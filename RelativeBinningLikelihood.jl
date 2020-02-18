#
# RelativeBinningLikelihood.jl
# Author: Nathaniel Leslie
# Relative Binning Paper: arXiv:1806.08792
#

# Packages not in PhenomDGradients

# Dependencies
include("PhenomD.jl")

# Constants
# ============================================================================
function RBLogLikelihoodGradient(LLOdata, LHOdata, LLOh0, LHOh0, LLOpsd, LHOpsd, θ, T, f_lo, f_hi, gmst, Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO)

    return gradient(MakeRBLogLikelihoodGradParams(LLOdata, LHOdata, LLOh0, LHOh0, LLOpsd, LHOpsd, T, f_lo, f_hi, gmst, Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO), θ)

end

function MakeRBLogLikelihoodGradParams(LLOdata, LHOdata, LLOh0, LHOh0, LLOpsd, LHOpsd, T, f_lo, f_hi, gmst, Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO)

    return function(θ) return RBLogLikelihood(LLOdata, LHOdata, LLOh0, LHOh0, LLOpsd, LHOpsd, θ, T, f_lo, f_hi, gmst, Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO) end

end

# ============================================================================
# Frequency array generator
function freqs(T, n_sample)
    fmax = n_sample/T
    return range(0.0, fmax/2.0, length=n_sample÷2+1)
end

# Find bin boundaries, bin indices, and number of bins, based on an improved version from the original paper, found in parameter_estimation_multidet.py (setup_data_binning_new)
# ✓
function GetBinParams(f, f_lo, f_hi, χ, ϵ)
    f_nt = range(f_lo, f_hi, length=10000) # nontrivial frequencies (where hp/hc are nonzero)

    γ = [-5/3, -2/3, 1, 5/3, 7/3]
    δα = 2*π*χ./abs.(f_lo.^γ-f_hi.^γ)
    dψ = zeros(length(f_nt))
    [dψ .+= sign(γ[i])*δα[i]*f_nt.^γ[i] for i=1:5]
    δΨ = dψ .- dψ[1]

    Nbin = convert(Int,δΨ[end]÷ϵ)
    Dϕ2f = LinearInterpolation(δΨ, f_nt)
    δΨgrid = range(δΨ[1], δΨ[end], length = Nbin+1)

    # Frequency grid points
    fbin = Dϕ2f(δΨgrid)

    # Find indices of frequencies in f
    fbin_ind = [findmin(abs.(f .- ff))[2] for ff in fbin]

    # Reset fbin to have those indices
    fbin = [f[ind] for ind in fbin_ind]

    # Remove repeated frequencies if any exist and set all params based on this list
    fbin = unique(fbin)
    Nbin = length(fbin) - 1
    fbin_ind = [findmin(abs.(f .- ff))[2] for ff in fbin]
    fm = 0.5*(fbin[2:end].+fbin[1:end-1])

    return Nbin, fbin, fbin_ind, fm
end

function ExactLogLikelihood(LLOdata, LHOdata, LLOpsd, LHOpsd, θ, T, n_sample, f_lo, f_hi, gmst)
    f = freqs(T, n_sample)

    Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance = θ

    hLLO, hLHO = PhenomDStrain(Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance, f, f_lo, f_hi, gmst)./T

    ZdhLLO, ZhhLLO = ComputeExactZOverlaps(LLOdata, hLLO, LLOpsd, T)
    ZdhLHO, ZhhLHO = ComputeExactZOverlaps(LHOdata, hLHO, LHOpsd, T)

    lnL_LLO = real(ZdhLLO) - 0.5*ZhhLLO
    lnL_LHO = real(ZdhLHO) - 0.5*ZhhLHO

    return lnL_LLO + lnL_LHO

end

function ComputeExactZOverlaps(d, h, psd, T)

    Zdh = 4*T*sum(d.*conj(h)./psd)
    Zhh = 4*T*sum(abs2.(h)./psd)

    return Zdh, Zhh
end

function SetupBinningData(LLOdata, LLOpsd, LHOdata, LHOpsd, fiducialθ, T, n_sample, f_lo, f_hi, χ, ϵ)

    f = freqs(T, n_sample)

    Mc_0, η_0, χ1_0, χ2_0, ra_0, dec_0, psi_0, iota_0, vphi_0, tc_0, distance_0 = fiducialθ
    LLOh0, LHOh0 = PhenomDStrain(Mc_0, η_0, χ1_0, χ2_0, ra_0, dec_0, psi_0, iota_0, vphi_0, tc_0, distance_0, f, f_lo, f_hi, gmst)./T

    Nbin, fbin, fbin_ind, fm = GetBinParams(f, f_lo, f_hi, χ, ϵ)

    sdatLLO = ComputeSummaryData(LLOdata, LLOh0, LLOpsd, T, f, Nbin, fbin, fbin_ind, fm)
    sdatLHO = ComputeSummaryData(LHOdata, LHOh0, LHOpsd, T, f, Nbin, fbin, fbin_ind, fm)

    return Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO, LLOh0, LHOh0

end

# The Log Likelihood is a constant offset from χ^2
function RBLogLikelihood(LLOdata, LHOdata, LLOh0, LHOh0, LLOpsd, LHOpsd, θ, T, f_lo, f_hi, gmst, Nbin, fbin, fbin_ind, fm, sdatLLO, sdatLHO)

    Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance = θ

    hLLO, hLHO = PhenomDStrain(Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance, fbin, f_lo, f_hi, gmst)./T

    rLLO = hLLO./LLOh0[fbin_ind]
    rLHO = hLHO./LHOh0[fbin_ind]

    r0LLO, r1LLO = ComputeBinCoefficients(fbin, rLLO)
    r0LHO, r1LHO = ComputeBinCoefficients(fbin, rLHO)

    ZdhLLO, ZhhLLO = ComputeZOverlaps(r0LLO, r1LLO, sdatLLO)
    ZdhLHO, ZhhLHO = ComputeZOverlaps(r0LHO, r1LHO, sdatLHO)

    lnL_LLO = real(ZdhLLO) - 0.5*ZhhLLO
    lnL_LHO = real(ZdhLHO) - 0.5*ZhhLHO

    return lnL_LLO + lnL_LHO
end

# The Log Likelihood is a constant offset from χ^2
# This version computes the fiducial waveform, bin parameters, and summary data on every call.
function RBLogLikelihood_slow(LLOdata, LHOdata, LLOpsd, LHOpsd, θ, fiducialθ, T, n_sample, f_lo, f_hi, gmst, χ, ϵ)

    f = freqs(T, n_sample)

    Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance = θ
    Mc_0, η_0, χ1_0, χ2_0, ra_0, dec_0, psi_0, iota_0, vphi_0, tc_0, distance_0 = fiducialθ

    Nbin, fbin, fbin_ind, fm = GetBinParams(f, f_lo, f_hi, χ, ϵ)

    hLLO, hLHO = PhenomDStrain(Mc, η, χ1, χ2, ra, dec, psi, iota, vphi, tc, distance, fbin, f_lo, f_hi, gmst)./T
    h0LLO, h0LHO = PhenomDStrain(Mc_0, η_0, χ1_0, χ2_0, ra_0, dec_0, psi_0, iota_0, vphi_0, tc_0, distance_0, f, f_lo, f_hi, gmst)./T

    rLLO = hLLO./h0LLO[fbin_ind]
    rLHO = hLHO./h0LHO[fbin_ind]

    r0LLO, r1LLO = ComputeBinCoefficients(fbin, rLLO)
    r0LHO, r1LHO = ComputeBinCoefficients(fbin, rLHO)

    sdatLLO = ComputeSummaryData(LLOdata, h0LLO, LLOpsd, T, f, Nbin, fbin, fbin_ind, fm)
    sdatLHO = ComputeSummaryData(LHOdata, h0LHO, LHOpsd, T, f, Nbin, fbin, fbin_ind, fm)

    ZdhLLO, ZhhLLO = ComputeZOverlaps(r0LLO, r1LLO, sdatLLO)
    ZdhLHO, ZhhLHO = ComputeZOverlaps(r0LHO, r1LHO, sdatLHO)

    lnL_LLO = real(ZdhLLO) - 0.5*ZhhLLO
    lnL_LHO = real(ZdhLHO) - 0.5*ZhhLHO

    return lnL_LLO + lnL_LHO
end

function ComputeZOverlaps(r0, r1, sdat)

    A₀, B₀, A₁, B₁ = sdat

    Zdh = sum(A₀.*conj(r0).+A₁.*conj(r1))
    Zhh = sum(B₀.*abs2.(r0)+2*B₁.*real(r0.*conj(r1)))

    return Zdh, Zhh
end

function ComputeSummaryData(d, h0, psd, T, f, Nbin, fbin, fbin_ind, fm)

    A₀ = 4*T*[sum(d[fbin_ind[i]:fbin_ind[i+1]-1].*conj(h0[fbin_ind[i]:fbin_ind[i+1]-1])./psd[fbin_ind[i]:fbin_ind[i+1]-1]) for i in range(1, stop=Nbin)]

    B₀ = 4*T*[sum(abs2.(h0[fbin_ind[i]:fbin_ind[i+1]-1])./psd[fbin_ind[i]:fbin_ind[i+1]-1]) for i in range(1, stop=Nbin)]

    A₁ = 4*T*[sum(d[fbin_ind[i]:fbin_ind[i+1]-1].*conj(h0[fbin_ind[i]:fbin_ind[i+1]-1])./psd[fbin_ind[i]:fbin_ind[i+1]-1].*(f[fbin_ind[i]:fbin_ind[i+1]-1].-fm[i])) for i in range(1, stop=Nbin)]

    B₁ = 4*T*[sum(abs2.(h0[fbin_ind[i]:fbin_ind[i+1]-1])./psd[fbin_ind[i]:fbin_ind[i+1]-1].*(f[fbin_ind[i]:fbin_ind[i+1]-1].-fm[i])) for i in range(1, stop=Nbin)]

    return A₀, B₀, A₁, B₁
end

function ComputeBinCoefficients(fbin, r)

    binwidths = fbin[2:end] - fbin[1:end-1]

    rplus = r[2:end]
    rminus = r[1:end-1]

    r0 = 0.5*(rplus + rminus)
    r1 = (rplus - rminus)./binwidths

    return r0, r1
end
