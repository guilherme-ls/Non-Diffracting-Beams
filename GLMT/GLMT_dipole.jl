using WGLMakie, ColorSchemes, ColorTypes, JLD2, SatelliteToolboxLegendre
using MakieThemes, SpecialFunctions, HypergeometricFunctions, AssociatedLegendrePolynomials

lambda = 1064 * 10^-9
c = 299792458
omega = 2 * pi * c / lambda
k = omega / c
K = 2 * k
Zmax = 22 * 10^-6
pmax = 1000
nmax = 47
Q = 0.9 * k
E0 = 1

z0 = 0
lmax = Int(round(Zmax * k / pi, RoundUp))
lmin = Int(round(-Zmax * k / pi, RoundDown))

alpha_var = 1 * 10^-3
a_mie = 1.0im * 2 / 3 * k^2 * alpha_var

tam = 50

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

@inline function F(z)
    term1 = 1 #exp(- 1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

Fterms = F.(pi .* range(lmin,lmax) ./ k)

function gTM(n, m)
    if abs(m) != 1
        return 0
    end

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
    if abs(m) != 1
        return 0
    end

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

@inline function G_factor()
    # not used, since m = 0 implies non-existant coefficients in the case of FWs
    # term1 = (1/3) * gTM(1, 0) * conj(gTM(2, 0))
    term2 = gTM(1, 1) * (conj(gTM(2, 1)) - 1im * conj(gTE(1, 1)))
    term3 = gTM(1, -1) * (conj(gTM(2, -1)) + 1im * conj(gTE(1, -1)))
    return term2 + term3
end

@inline function Fz()
    term1 = 3 * lambda^2 / (2 * pi)
    term2 = real(a_mie * G_factor())
    return term1 * term2
end

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Fzvalues = []
global count = 0
for z in zvalues
    global z0 = z
    push!(Fzvalues, Fz())
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
    xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (N)", aspect=1.5)

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
