using JLD2
using SpecialFunctions, HypergeometricFunctions

pmax = 1000
nmax = 50

function Qn(n,m)
    if m >= n
        return 0
    end
    return trunc(Int, (n-m)/2)
end

function a(n, m, q) # uses log to prevent overflow with big factorials
    log_num = loggamma(n + m + 1)
    log_denom = loggamma(m + q + 1) +
                loggamma(n - m - 2q + 1) +
                loggamma(q + 1) +
                (m + 2q) * log(2)

    return (-1.0)^q * exp(log_num - log_denom)
end

function B1(n,p)
    term1 = n * (n + 1) / 2
    term2 = 0
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = (n - 2q + 1)/2
            term2 += - a(n,0,q) * beta(mu, nu) * pFq((nu,), (0.5, nu+mu), - pi^2 * p^2 / 4)
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = (n - 2q + 1)/2
            term2 += a(n,0,q) * 1im * pi * p * beta(mu, nu + 0.5) * pFq((nu + 0.5,), (1.5, nu + mu + 0.5), - pi^2 * p^2 / 4)
        end
    end
    return term1 * term2
end

function B2(n, p)
    return B1(n,p) / (n * (n+1))
end

function B3A(n,p)
    term1 = - 0.5
    term2 = 0
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,1))
            mu = q + 2
            nu = (n - q)/2
            term2 += - a(n,1,q) * beta(mu, nu) * pFq((nu,), (0.5, nu+mu), - pi^2 * p^2 / 4)
        end
    else
        for q in range(0, Qn(n,1))
            mu = q + 2
            nu = (n - q)/2
            term2 += a(n,1,q) * 1im * pi * p * beta(mu, nu + 0.5) * pFq((nu + 0.5,), (1.5, nu + mu + 0.5), - pi^2 * p^2 / 4)
        end
    end
    return term1 * term2
end

function B3B(n,p)
    term1 = 0.5
    term2 = 0
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = (n - 2q + 2)/2
            term2 += - a(n,0,q) * beta(mu, nu) * pFq((nu,), (0.5, nu+mu), - pi^2 * p^2 / 4)
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = (n - 2q + 2)/2
            term2 += a(n,0,q) * 1im * pi * p * beta(mu, nu + 0.5) * pFq((nu + 0.5,), (1.5, nu + mu + 0.5), - pi^2 * p^2 / 4)
        end
    end
    return term1 * term2
end

function B3(n,p)
    return B3A(n,p) + B3B(n,p)
end

function B4(n, p)
    return B3(n,p) / (n * (n+1))
end


function save_B_values(filename::String, nmax::Int, pmax::Int)
    B1_vals = Dict{Tuple{Int, Int}, ComplexF64}()
    B2_vals = Dict{Tuple{Int, Int}, ComplexF64}()
    B3_vals = Dict{Tuple{Int, Int}, ComplexF64}()
    B4_vals = Dict{Tuple{Int, Int}, ComplexF64}()

    for n in 1:nmax
        for p in -pmax:pmax
            key = (n, p)
            B1_vals[key] = B1(n, p)
            B2_vals[key] = B2(n, p)
            B3_vals[key] = B3(n, p)
            B4_vals[key] = B4(n, p)
        end
    end

    @save filename B1_vals B2_vals B3_vals B4_vals
end

save_B_values("./GLMT/Bvalues.jld2", nmax, pmax)
