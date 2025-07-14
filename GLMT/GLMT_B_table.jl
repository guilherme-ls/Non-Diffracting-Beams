using JLD2, Base.Threads
using SpecialFunctions, HypergeometricFunctions, Arblib

pmax = 500
nmax = 80

setprecision(256)
Arblib.setprecision(256)

function Qn(n,m)
    if m >= n
        return 0
    end
    return trunc(Int, (n-m)/2)
end

function a(n, m, q) # uses log to prevent overflow with big factorials
    log_num = loggamma(BigFloat(n + m + 1))
    log_denom = loggamma(BigFloat(m + q + 1)) + loggamma(BigFloat(n - m - 2q + 1)) + loggamma(BigFloat(q + 1)) + (m + 2q) * log(BigFloat(2))

    return (-1.0)^q * exp(log_num - log_denom)
end

function hypergeometric_1F2(a::BigFloat, b1::BigFloat, b2::BigFloat, z::BigFloat)
    return pFq((a,), (b1,b2),z)
end

Aterms = [a(n,m,q) for n in 1:nmax, m in 0:1, q in 0:(nmax/2)]
println("Calculated Aterms")

Bterms1 = [beta(BigFloat(mu), BigFloat(nu/2)) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 2)]
Bterms2 = [beta(BigFloat(mu), BigFloat((nu + 1)/2)) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 2)]
println("Calculated Bterms")

pFqterms1 = Array{Complex{BigFloat}}(undef, Int(round(nmax/2 + 2, RoundDown)), nmax + 2, pmax * 2 + 1)
pFqterms2 = Array{Complex{BigFloat}}(undef, Int(round(nmax/2 + 2, RoundDown)), nmax + 2, pmax * 2 + 1)
start_time = time()
for p in 0:pmax
    param = - big(pi)^2 * p^2 / 4
    Threads.@threads for mu in 1:Int(round(nmax/2 + 2, RoundDown))
        for nu in 1:(nmax + 2)
            #pFqterms1[mu, nu, p + pmax + 1] = pFq((BigFloat(nu/2),), (BigFloat(0.5), BigFloat(nu/2+mu)), param)
            #pFqterms2[mu, nu, p + pmax + 1] = pFq((BigFloat((nu + 1)/2),), (BigFloat(1.5), BigFloat((nu + 1)/2 + mu)), param)
            pFqterms1[mu, nu, p + pmax + 1] = pFqterms1[mu, nu, -p + pmax + 1] = hypergeometric_1F2(BigFloat(nu/2), BigFloat(0.5), BigFloat(nu/2+mu), param)
            pFqterms2[mu, nu, p + pmax + 1] = pFqterms2[mu, nu, -p + pmax + 1] = hypergeometric_1F2(BigFloat((nu + 1)/2), BigFloat(1.5), BigFloat((nu + 1)/2 + mu), param)
        end
    end
    if p%10 == 0
        println("Calculou ", 2 * p, " termos de pFq. Tempo decorrido: ", time() - start_time)
    end
end
println("Calculated pFqterms")

function B1(n,p)
    term1 = BigFloat(n * (n + 1) / 2)
    term2 = BigFloat(0)
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += - Aterms[n,1,q+1] * Bterms1[mu, nu] * pFqterms1[mu, nu, p + pmax + 1]
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += Aterms[n,1,q+1] * 1im * big(pi) * p * Bterms2[mu, nu] * pFqterms2[mu, nu, p + pmax + 1]
        end
    end
    return term1 * term2
end

function B2(n, p)
    return B1(n,p) / BigFloat((n * (n+1)))
end

function B3A(n,p)
    term1 = BigFloat(- 0.5)
    term2 = BigFloat(0)
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,1))
            mu = q + 2
            nu = n - q
            term2 += - Aterms[n,2,q+1] * Bterms1[mu, nu] * pFqterms1[mu, nu, p + pmax + 1]
        end
    else
        for q in range(0, Qn(n,1))
            mu = q + 2
            nu = n - q
            term2 += Aterms[n,2,q+1] * 1im * big(pi) * p * Bterms2[mu, nu] * pFqterms2[mu, nu, p + pmax + 1]
        end
    end
    return term1 * term2
end

function B3B(n,p)
    term1 = BigFloat(0.5)
    term2 = BigFloat(0)
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 2
            term2 += - Aterms[n,1,q+1] * Bterms1[mu, nu] * pFqterms1[mu, nu, p + pmax + 1]
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 2
            term2 += Aterms[n,1,q+1] * 1im * big(pi) * p * Bterms2[mu, nu] * pFqterms2[mu, nu, p + pmax + 1]
        end
    end
    return term1 * term2
end

function B3(n,p)
    return B3A(n,p) + B3B(n,p)
end

function B4(n, p)
    return B3(n,p) / BigFloat(n * (n+1))
end


function save_B_values(filename::String)
    B1_vals = Array{Complex{BigFloat}}(undef, nmax, pmax * 2 + 1)
    B2_vals = Array{Complex{BigFloat}}(undef, nmax, pmax * 2 + 1)
    B3_vals = Array{Complex{BigFloat}}(undef, nmax, pmax * 2 + 1)
    B4_vals = Array{Complex{BigFloat}}(undef, nmax, pmax * 2 + 1)

    for n in 1:nmax
        for p in -pmax:pmax
            B1_vals[n, p + pmax + 1] = B1(n, p)
            B2_vals[n, p + pmax + 1] = B1_vals[n, p + pmax + 1] / BigFloat(n * (n+1))
            B3_vals[n, p + pmax + 1] = B3(n, p)
            B4_vals[n, p + pmax + 1] = B3_vals[n, p + pmax + 1] / BigFloat(n * (n+1))
        end
    end

    @save filename B1_vals B2_vals B3_vals B4_vals
end

save_B_values("./GLMT/Bvalues.jld2")

println("B1(1, 1): ", B1(1,1))
println("B1(2, 1): ", B1(2,1))
println("B1(15, 1): ", B1(15,1))
println("B1(20, 1): ", B1(20,1))
println("B1(25, 1): ", B1(25,1))