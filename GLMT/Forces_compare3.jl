using GLMakie, ColorSchemes, ColorTypes, JLD2, FiniteDifferences
using MakieThemes, SpecialFunctions, AssociatedLegendrePolynomials, Arblib, Integrals

setprecision(256)
Arblib.setprecision(256)

lambda = BigFloat("1064e-9")
mu = 4 * big(pi) * BigFloat("1e-7")
epsilon = BigFloat(885418782.0) * BigFloat("1e-20")
c = BigFloat(1 / sqrt(mu * epsilon))
k = 2 * big(pi) / lambda
omega = k * c
K = 2 * k
z0 = 0
Q = 0.8 * k
E0 = 1
eta_var = sqrt(mu / epsilon)
mult_factor = BigFloat("1e+24")

M = 1.2
d = 10^-9
alpha = big(pi) * d / lambda
betah = M * alpha
musp = mu

pmax = 500
nmax = 50

Zmax = 40 * BigFloat(10^-6)
lmax = Int(round(Zmax * k / pi, RoundDown))
lmin = Int(round(-Zmax * k / pi, RoundUp))

tam = 100

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

function F(z)
    term1 = exp(Complex{BigFloat}(0, 1) * Q * z)
    term2 = exp(- 8 * (z/Zmax)^8)
    return term1 * term2
end

function f(z, p)
    return F(z)
end

function g(z)
    prob = IntegralProblem(f, (-1.2 * Zmax, z))
    return solve(prob, HCubatureJL(), reltol = 1e-6, abstol = 1e-6).u
end

Fterms = F.(- big(pi) .* range(lmin,lmax) ./ k)

function custom_sinint(x)
    x_arb = Arb(x)
    si = sinint(x_arb)
    return BigFloat(si)
end

sinint_table = []
function gTM(n, m)
    term1 = - (1.0im)^m / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = BigFloat(0)
    for p in range(-pmax,pmax)
        term3 = BigFloat(0)
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint_table[l - lmin + p + pmax + 1]
        end
        term2 += (B1_data[n, p + pmax + 1] * (m + 1) + B2_data[n, p + pmax + 1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function gTE(n, m)
    term1 = (1.0im)^(m + 1) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint_table[l - lmin + p + pmax + 1]
        end
        term2 += (B3_data[n, p + pmax + 1] * (m + 1) - B4_data[n, p + pmax + 1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function custom_besselj(n, x)
    n_acb = Acb(n + 0.5)
    x_acb = Acb(x)
    bess = besselj(n_acb, x_acb)
    result = sqrt(big(pi) / (2 * x)) * Complex{BigFloat}(bess)
    return result
end

function custom_hankel2(n, x)
    n_acb = Acb(n + 0.5)
    x_acb = Acb(x)
    hank = besselh(n_acb, 2, x_acb)
    result = sqrt(big(pi) / (2 * x)) * Complex{BigFloat}(hank)
    return result
end

function ricatti_besselj(n, x)
    return x * custom_besselj(n, x)
end

function dricatti_besselj(n, x)
    #return forward_fdm(5,1)(x1 -> ricatti_besselj(n, x1), x)
    return - n * custom_besselj(n, x) + x * custom_besselj(n-1, x)
    #return (n + 1) * custom_besselj(n, x) - x * custom_besselj(n-1, x)
end

function ricatti_besselh(n, x)
    return x * custom_hankel2(n, x)
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
    term1 = mu * M * ricatti_besselj(n, alpha) * dricatti_besselj(n, betah) - 
        musp * ricatti_besselj(n, betah) * dricatti_besselj(n, alpha)
    term2 = mu * M * ricatti_besselh(n, alpha) * dricatti_besselj(n, betah) - 
        musp * ricatti_besselj(n, betah) * dricatti_besselh(n, alpha)
    return term1 / term2
end

nterms = range(1,nmax)
anterms = an.(nterms)
bnterms = bn.(nterms)

function Cpr_z()
    C = BigFloat(0)

    Threads.@threads for n in 1:5
        for p in (-1,1) # the only valid ones
            abs_p = abs(p)

            # Factorial prefactor
            B = exp(loggamma(BigFloat(n + abs_p + 1)) - loggamma(BigFloat(n - abs_p + 1)))
            A = B * (n + 1 + abs_p)

            # precalculates relevant LMT coefficients
            const_a = anterms[n]
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
            T1_real = real((const_a + conj(const_a1) - 2 * const_a * conj(const_a1)) * gtma * conj(gtma1) +
                (const_b + conj(const_b1) - 2 * const_b * conj(const_b1)) * gtea * conj(gtea1))

            # Second Term
            T2_factor = p * BigFloat(2n + 1) / BigFloat(n^2 * (n + 1)^2) * B
            T2_real = real(1im * (2 * const_a * conj(const_b) - const_a - conj(const_b)) * gtma * conj(gtea))

            C += T1_factor * T1_real + T2_factor * T2_real
        end
    end

    return C * lambda^2 / big(pi) * mult_factor
end

function Ex(r, z)
    values = 0
    for n in range(lmin,lmax)
        zf = - big(pi) * n / k
        values += F(zf) * sinc.(1/big(pi) * sqrt((k * r)^2 + (z * k + big(pi) * n)^2))
    end
    return E0 * values
end

function Ez(r, z, phi)
    values = 0
    for n in range(lmin,lmax)
        zf = - big(pi) * n / k
        values += g(zf) * k^2 * r * (cos.(sqrt.((k * r)^2 + (z * k + pi * n)^2))/((k * r)^2 + (z * k + pi * n)^2) - sin.(sqrt.((k * r)^2 .+ (z * k .+ pi .* n).^2))./((k * r)^2 .+ (z * k .+ pi .* n).^2).^1.5)
    end
    return -E0 * cos(phi) * values
end

function dsinc_analytic(r, z, n)
    factor = k^2 * r^2 + (k * z + big(pi) * n)^2
    sqrtA = sqrt(factor)
    
    if factor == 0.0
        return 0.0
    end

    term1 = cos(sqrtA) / sqrtA
    term2 = sin(sqrtA) / factor
    prefactor = k * (k * z + big(pi) * n)

    return prefactor * (term1 - term2) / sqrtA
end

function dEx_conj(r, z)
    values = 0
    for n in range(lmin,lmax)
        zf = - big(pi) * n / k
        values += F(zf) * dsinc_analytic(r, z, n)
    end
    return conj(E0 * values)
end

function dEz_conj(r, z, phi)
    return conj(central_fdm(3,1)(z1 -> Ez(r, z1, phi), z))
end

a_mie = an(1)
alpha_var = - 1im * 3 / (2 * k^3) * a_mie

function Fz_dipole(r, z, phi)
    term1 = 2 * big(pi) / eta_var
    term2 = alpha_var * (Ex(r,z) * dEx_conj(r,z) + Ez(r,z,phi) * dEz_conj(r,z,phi))
    return term1 * real(term2) * mult_factor
end

@inline function G_factor()
    term1 = (1/3) * gTM(1, 0) * conj(gTM(2, 0))
    term2 = gTM(1, 1) * (conj(gTM(2, 1)) - 1im * conj(gTE(1, 1)))
    term3 = gTM(1, -1) * (conj(gTM(2, -1)) + 1im * conj(gTE(1, -1)))
    return term1 + term2 + term3
end

@inline function Fz_GLMT()
    term1 = 3 * lambda^2 / (2 * big(pi))
    term2 = real(a_mie * G_factor())
    return term1 * term2 * mult_factor
end

#############################################################
# ALL CALCULATIONS ARE DONE HERE
#############################################################

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Fzvalues_dipole = Array{BigFloat}(undef, tam)
Threads.@threads for count in 1:tam
    global Fzvalues_dipole[count] =  Fz_dipole(0, -zvalues[count], 0)
    if count % 10 == 0
        println(count, " pontos calculados")
    end
end

Fzvalues_GLMT = Array{BigFloat}(undef, tam)
Cvalues = []
global count = 0
for z in zvalues
    global z0 = z
    global sinint_table = [custom_sinint(big(pi) * term + k * z0) for term in (-pmax+lmin):(pmax+lmax)]
    push!(Cvalues, Cpr_z())
    global Fzvalues_GLMT[count + 1] = Fz_GLMT()
    if count % 10 == 0
        println(count, " pontos calculados")
    end
    global count += 1
end
z0 = 0

function FW(r, z)
    sum_mackinnon = 0 + 0im
    for n in range(lmin, lmax)
        sum_mackinnon += F(- n * pi / k) * sinc(1/pi * sqrt((k * r)^2 + (z * k + pi * n)^2))
    end
    return sum_mackinnon
end

function S(x)
    return 1/K * reduce(+, Fterms .* exp.(1im * 2 * pi * range(lmin,lmax) * x / K))
end

intensityFW = [abs2(FW(0, z)) for z in zvalues]
spectra = [abs2(S(x)) for x in range(-k,k,tam)]
println("FW calculated")

# plots
fig = Figure(figure_padding=40, size= 3 .*(800, 200), fontsize = 25)
ax1 = Axis(fig[1, 1], title = "(a)", xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.3)
ax2 = Axis(fig[1, 2], title = "(b)", yticklabelsvisible=false, yticksvisible=false,
    xgridvisible = false, ygridvisible = false, xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (pm^2)", aspect=1.3)
ax3 = Axis(fig[1, 3], title = "(c)", yticklabelsvisible=false, yticksvisible=false,
    xgridvisible = false, ygridvisible = false, ylabelvisible=false, xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (pm^2)", aspect=1.3)
ax4 = Axis(fig[1, 4], title = "(d)", yticklabelsvisible=false, yticksvisible=false,
    xgridvisible = false, ygridvisible = false, ylabelvisible=false, xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (pm^2)", aspect=1.3)

ax_inset = Axis(fig[1, 1],
    width=Relative(0.2),
    height=Relative(0.2),
    halign=0.05, valign=0.85,
    xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
    ylabelvisible=false, xticksvisible=false)

Fplot = abs2.(F.(zvalues))
zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="black", label=L"|F(z)|^2")
lines!(ax1, zvalues, intensityFW, color="blue", label=L"|\Psi(0,z)|^2", linestyle=:dash)
lines!(ax_inset, zvalues, spectra, color="black")
lines!(ax2, zvalues, Fzvalues_dipole, color="black")
lines!(ax3, zvalues, Cvalues, color="black")
lines!(ax4, zvalues, Fzvalues_GLMT, color="black")

axislegend(ax1, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("Forces_compare3.png", fig)
