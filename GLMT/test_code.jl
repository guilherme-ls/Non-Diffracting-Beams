using JLD2, FiniteDifferences, Arblib
using SpecialFunctions, HypergeometricFunctions, AssociatedLegendrePolynomials

setprecision(256)
Arblib.setprecision(256)

lambda = BigFloat("1064e-9")
mu = 4 * big(pi) * BigFloat("1e-7")
epsilon = 885418782 * BigFloat("1e-20")
c = BigFloat(1 / sqrt(mu * epsilon))
k = 2 * big(pi) / lambda
omega = k * c
K = 2 * k
Zmax = 22 * BigFloat(10^-6)
z0 = 0
Q = 0.9 * k
E0 = 1

lmax = Int(round(Zmax * k / pi, RoundUp))
lmin = Int(round(-Zmax * k / pi, RoundDown))

pmax = 500
nmax = 47

function F(z)
    term1 = exp(- 1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

Fterms = F.(big(pi) .* range(lmin,lmax) ./ k)

function Qn(n,m)
    if m >= n
        return 0
    end
    return Int(round((n-m)/2, RoundDown))
end

function a(n, m, q) # uses log to prevent overflow with big factorials
    log_num = loggamma(BigFloat(n + m + 1))
    log_denom = loggamma(BigFloat(m + q + 1)) + loggamma(BigFloat(n - m - 2q + 1)) + loggamma(BigFloat(q + 1)) + (m + 2q) * log(BigFloat(2))

    return (-1.0)^q * exp(log_num - log_denom)
end

function hypergeometric_1F2(a::BigFloat, b1::BigFloat, b2::BigFloat, z::BigFloat)
    if z > -6000
        a_acb = AcbVector([a])
        b_acb = AcbVector([b1,b2])
        z_acb = Acb(z)
        result = Acb()
        Arblib.hypgeom_pfq!(result, a_acb, 1, b_acb, 2, z_acb, 0)
        return Complex{BigFloat}(result)
    else
        return pFq((a,), (b1,b2),z)
    end
end

Aterms = [a(n,m,q) for n in 1:nmax, m in 0:1, q in 0:(nmax/2)]
println("Calculated Aterms")

Bterms1 = [beta(BigFloat(mu), BigFloat(nu/2)) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 2)]
Bterms2 = [beta(BigFloat(mu), BigFloat((nu + 1)/2)) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 2)]
println("Calculated Bterms")

#pFqterms1 = [pFq((BigFloat(nu/2),), (BigFloat(0.5), BigFloat(nu/2+mu)), - big(pi)^2 * p^2 / 4) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 1), p in -pmax:pmax]
println("Calculated pFqterms1")

#pFqterms2 = [pFq((BigFloat((nu + 1)/2),), (BigFloat(1.5), BigFloat((nu + 1)/2 + mu)), - big(pi)^2 * p^2 / 4) for mu in 1:(nmax/2 + 2), nu in 1:(nmax + 1), p in -pmax:pmax]
println("Calculated pFqterms2")

function B1(n,p)
    term1 = BigFloat(n * (n + 1) / 2)
    term2 = BigFloat(0)
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += - Aterms[n,1,q+1] * Bterms1[mu, nu] * hypergeometric_1F2(BigFloat(nu/2), BigFloat(0.5), BigFloat(nu/2+mu), - big(pi)^2 * p^2 / 4)
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += Aterms[n,1,q+1] * 1im * big(pi) * p * Bterms2[mu, nu] * hypergeometric_1F2(BigFloat((nu + 1)/2), BigFloat(1.5), BigFloat((nu + 1)/2 + mu), - big(pi)^2 * p^2 / 4)
        end
    end
    return term1 * term2
end

function B1_alt(n,p)
    term1 = BigFloat(n * (n + 1) / 2)
    term2 = BigFloat(0)
    if n%2 == 0 # checks if n is even or odd
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += - Aterms[n,1,q+1] * Bterms1[mu, nu] * pFq((BigFloat(nu/2),), (BigFloat(0.5), BigFloat(nu/2+mu)), - big(pi)^2 * p^2 / 4)
        end
    else
        for q in range(0, Qn(n,0))
            mu = q + 1
            nu = n - 2q + 1
            term2 += Aterms[n,1,q+1] * 1im * big(pi) * p * Bterms2[mu, nu] * pFq((BigFloat((nu + 1)/2),), (BigFloat(1.5), BigFloat((nu + 1)/2 + mu)), - big(pi)^2 * p^2 / 4)
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

#B1_data = [B1(n,p) for n in 1:nmax, p in -pmax:pmax]
#B2_data = [B1_data[n,p + pmax +1] / BigFloat((n * (n+1))) for n in 1:nmax, p in -pmax:pmax]
println("Calculated B1 and B2 terms")
#B3_data = [B3(n,p) for n in 1:nmax, p in -pmax:pmax]
#B4_data = [B3_data[n,p + pmax +1] / BigFloat((n * (n+1))) for n in 1:nmax, p in -pmax:pmax]
println("Calculated B3 and B4 terms")

function gTM(n, m)
    if abs(m) != 1
        return 0
    end

    term1 = - (1.0im)^m / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = BigFloat(0)
    for p in range(-pmax,pmax)
        term3 = BigFloat(0)
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * BigFloat(sinint(Float64(pi * (l + p) + k * z0)))
        end
        term2 += (B1_data[n, p+pmax+1] * (m + 1) + B2_data[n, p+pmax+1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function gTE(n, m)
    if abs(m) != 1
        return 0
    end

    term1 = (1.0im)^(m + 1) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint(pi * (l + p) + k * z0)
        end
        term2 += (B3_data[n, p+pmax+1] * (m + 1) - B4_data[n, p+pmax+1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

point = -40

println("B1(1, 1): ", B1(1,point))
println("B1(3, 1): ", B1(3,point))
println("B1(15, 1): ", B1(15,point))
println("B1(20, 1): ", B1(20,point))
println("B1(25, 1): ", B1(25,point))

println("B1(1, 1) esperado: ", B1_alt(1,point))
println("B1(3, 1) esperado: ", B1_alt(3,point))
println("B1(15, 1) esperado: ", B1_alt(15,point))
println("B1(20, 1) esperado: ", B1_alt(20,point))
println("B1(25, 1) esperado: ", B1_alt(25,point))

#println("gTM n=1: ", gTM(1,1))
#println("gTM n=2: ", gTM(2,1))
#println("gTM n=5: ", gTM(5,1))
#println("gTM n=10: ", gTM(10,1))
#println("gTM n=25: ", gTM(25,1))