using GLMakie, ColorSchemes, ColorTypes, JLD2, Base.Threads, FiniteDifferences
using MakieThemes, SpecialFunctions, HypergeometricFunctions, AssociatedLegendrePolynomials, Arblib

setprecision(256)
Arblib.setprecision(256)

lambda = BigFloat("1064e-9")
mu = 4 * big(pi) * BigFloat("1e-7")
epsilon = 885418782 * BigFloat("1e-20")
c = BigFloat(1 / sqrt(mu * epsilon))
k = 2 * big(pi) / lambda
omega = k * c
K = 2 * k
Zmax = 30 * BigFloat(10^-6)
z0 = 0
Q = 0.9 * k
E0 = 1
mult_factor = BigFloat("1e+24")

pmax = 500
nmax = 47

z0 = 0
lmax = Int(round(Zmax * k / pi, RoundUp))
lmin = Int(round(-Zmax * k / pi, RoundDown))

tam = 500

M = 1.2
d = 10^-9
alpha = big(pi) * d / lambda
betah = M * alpha
musp = mu

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
a_mie = an(1)
println(a_mie)

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

function F(z)
    term1 = exp(1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

Fterms = F.(- big(pi) .* range(lmin,lmax) ./ k)

function custom_sinint(x)
    x_arb = Arb(x)
    si = sinint(x_arb)
    return BigFloat(si)
end

prev_z0 = Inf
sinint_table = Array{BigFloat}(undef, 2 * pmax - lmin + lmax + 1)
function gTM(n, m)
    if abs(m) != 1
        return 0
    end

    if prev_z0 != z0
        global prev_z0 = z0
        Threads.@threads for term in (-pmax+lmin):(pmax+lmax)
            global sinint_table[term - lmin + pmax + 1] = custom_sinint(big(pi) * term + k * z0)
        end
    end

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
    if abs(m) != 1
        return 0
    end

    if prev_z0 != z0
        global prev_z0 = z0
        Threads.@threads for term in (-pmax+lmin):(pmax+lmax)
            global sinint_table[term - lmin + pmax + 1] = custom_sinint(big(pi) * term + k * z0)
        end
    end

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

@inline function G_factor()
    term1 = (1/3) * gTM(1, 0) * conj(gTM(2, 0))
    term2 = gTM(1, 1) * (conj(gTM(2, 1)) - 1im * conj(gTE(1, 1)))
    term3 = gTM(1, -1) * (conj(gTM(2, -1)) + 1im * conj(gTE(1, -1)))
    return term1 + term2 + term3
end

@inline function Fz()
    term1 = 3 * lambda^2 / (2 * big(pi))
    term2 = real(a_mie * G_factor())
    return term1 * term2 * mult_factor
end

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Fzvalues = Array{BigFloat}(undef, tam)
global count = 0
for z in zvalues
    global z0 = z
    global Fzvalues[count + 1] = Fz()
    if count % 100 == 0
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
    xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (pm^2)", aspect=1.5)

Fplot = abs2.(F.(zvalues))

zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="blue")
lines!(ax2, zvalues, Fzvalues, color="blue")
#lines!(ax2, zvalues, Fzvalues, color="black")

#axislegend(ax1, position=:rt)
#axislegend(ax2, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("GLMT_dipole.png", fig)
