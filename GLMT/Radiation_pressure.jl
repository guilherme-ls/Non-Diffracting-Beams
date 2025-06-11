using WGLMakie, ColorSchemes, ColorTypes, JLD2, FiniteDifferences
using MakieThemes, SpecialFunctions, Interpolations, Integrals, HypergeometricFunctions, AssociatedLegendrePolynomials

lambda = 1064 * 10^-9
c = 299792458
omega = 2 * pi * c / lambda
k = omega / c
K = 2 * k
Zmax = 22 * 10^-6
pmax = 500
nmax = 47
global z0 = 0
Q = 0.9 * k
E0 = 1

mu = 1.256 * 10^-6
musp = mu
M = 1.59
d = 10^-6
alpha = pi * d / lambda
betah = M * alpha

lmax = Int(round(Zmax * k / pi, RoundUp))
lmin = Int(round(-Zmax * k / pi, RoundDown))

tam = 24

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

@inline function F(z)
    term1 = 1 # exp(- 1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

Fterms = F.(pi .* range(lmin,lmax) ./ k)

function gTM(n, m)
    term1 = - (1.0im)^m / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(n-m+1) - loggamma(n + abs(m) + 1)) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint(pi * (l + p) + k * z0)
        end
        term2 += (B1_data[(n, p)] * (m + 1) + B2_data[(n, p)] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function gTE(n, m)
    term1 = (1.0im)^(m + 1) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(n-m+1) - loggamma(n + abs(m) + 1)) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint(pi * (l + p) + k * z0)
        end
        term2 += (B3_data[(n, p)] * (m + 1) - B4_data[(n, p)] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function ricatti_besselj(n, x)
    return x * sphericalbesselj(n, x)
end

function dricatti_besselj(n, x)
    return forward_fdm(5,1)(x1 -> ricatti_besselj(n, x1), x)
end

function sphericalhankel(n, x)
    return sqrt(pi / (2 * x)) * hankelh2(n + 0.5, x)
end

function ricatti_besselh(n, x)
    return x * sphericalhankel(n, x)
end

function dricatti_besselh(n, x)
    return forward_fdm(5,1)(x1 -> ricatti_besselh(n, x1), x)
end

function an(n)
    term1 = musp * ricatti_besselj(n, alpha) * dricatti_besselj(n, betah) - 
        mu * M * ricatti_besselj(n, betah) * dricatti_besselj(n, alpha)
    term2 = musp * ricatti_besselh(n, alpha) * dricatti_besselj(n, betah) - 
        mu * M * ricatti_besselj(n, betah) * dricatti_besselh(n, alpha)
    return term1 / term2
end

function bn(n)
    term1 = mu * ricatti_besselj(n, alpha) * dricatti_besselj(n, betah) - 
        musp * M * ricatti_besselj(n, betah) * dricatti_besselj(n, alpha)
    term2 = mu * ricatti_besselh(n, alpha) * dricatti_besselj(n, betah)- 
        musp * M * ricatti_besselj(n, betah) * dricatti_besselh(n, alpha)
    return term1 / term2
end

nterms = range(1,nmax)
anterms = an.(nterms)
bnterms = bn.(nterms)

function Cpr_z()
    C = 0.0

    for n in 1:(nmax-1)
        for p in (-1,1) # the only valid ones
            abs_p = abs(p)

            # Factorial prefactor
            B = exp(loggamma(n + abs_p + 1) - loggamma(n - abs_p + 1))
            A = B * (n + 1 + abs_p)

            # precalculates relevant LMT coefficients
            const_a = anterms[n+1]
            const_a1 = anterms[n+1]
            const_b = bnterms[n]
            const_b1 = bnterms[n+1]

            # precalculates relevant g terms
            gtma = gTM(n, p)
            gtma1 = gTM(n+1, p)
            gtea = gTE(n, p)
            gtea1 = gTE(n+1, p)

            # First Term
            T1_factor = 1 / (n + 1)^2 * A
            T1_real = real((const_a + conj(const_a1) - 2const_a * conj(const_a1)) * gtma * conj(gtma1) +
                (const_b + conj(const_b1) - 2const_b * conj(const_b1)) * gtea * conj(gtea1))

            # Second Term
            T2_factor = p * (2n + 1) / (n^2 * (n + 1)^2) * B
            T2_real = real(1im * (2const_a * conj(const_b) - const_a - conj(const_b)) * gtma * conj(gtea))

            C += T1_factor * T1_real + T2_factor * T2_real
        end
    end

    return lambda^2 / pi * C
end

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Cvalues = []
global count = 0
for z in zvalues
    global z0 = z
    push!(Cvalues, Cpr_z())
    if count % 10 == 0
        println(count, " pontos calculados")
    end
    global count += 1
end
z0 = 0

# plots
fig = Figure(figure_padding=40, size= 3 .*(500, 200), fontsize = 25)
ax1 = Axis(fig[1, 1], title = "(a)",
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.5)
ax2 = Axis(fig[1, 2], title = "(b)",
    xlabel = L"z_0 \, (\mu m)", ylabel = L"C_{pr,z} \, (m^2)", aspect=1.5)

Fplot = abs2.(F.(zvalues))

zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="blue")
lines!(ax2, zvalues, Cvalues, color="black")

#axislegend(ax1, position=:rt)
#axislegend(ax2, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("Radiation_pressure.png", fig)
