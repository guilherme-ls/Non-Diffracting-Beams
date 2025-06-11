using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions, Interpolations, Integrals, FiniteDifferences

lambda = 1064 * 10^-9
c = 299792458
omega = 2 * pi * c / lambda
k = omega / c
K = 2 * k
Zmax = 22 * 10^-6
#N = Int(round(Z * omega / (pi * c), RoundDown))
Q = 0.9 * k
E0 = 1

alpha_var = 1 * 10^-3
eta_var = 337
z0 = 0
nmax = Int(round((Zmax + z0) * k / pi, RoundUp))
nmin = Int(round((-Zmax + z0) * k / pi, RoundDown))

tam = 1000

@inline function F(z)
    term1 = 1 # exp(- 1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

function g(z)
    values = LinRange(-10 * Z, -2 * pi * n /K, 100000)
    soma = reduce(+, F.(values) * (values[2] - values[1]))
    return soma
end

function Ex(r, z)
    values = 0
    for n in range(nmin,nmax)
        zf = pi * n / k
        values += F(zf) * sinc.(1/pi * sqrt((k * r)^2 + (z * k + pi * n)^2))
    end
    return E0 * values
end

function Ez(r, z, phi)
    if r == 0
        return 0
    end
    values = 0
    for n in range(nmin,nmax)
        zf = pi * n / k
        values += g(zf) * k^2 * r * (cos.(sqrt.((k * r)^2 + (z * k + pi * n)^2))/((k * r)^2 + (z * k + pi * n)^2) - sin.(sqrt.((k * r)^2 .+ (z * k .+ pi .* n).^2))./((k * r)^2 .+ (z * k .+ pi .* n).^2).^1.5)
    end
    return -E0 * cos(phi) * values
end

function dsinc_analytic(r, z, n)
    factor = k^2 * r^2 + (k * z + pi * n)^2
    sqrtA = sqrt(factor)
    
    if factor == 0.0
        return 0.0
    end

    term1 = cos(sqrtA) / sqrtA
    term2 = sin(sqrtA) / factor
    prefactor = k * (k * z + pi * n)

    return prefactor * (term1 - term2) / sqrtA
end

function dEx_conj(r, z)
    values = 0
    for n in range(nmin,nmax)
        zf = - pi * n / k
        values += F(zf) * dsinc_analytic(r, z, n)
    end
    return conj(E0 * values)
end

function dEz_conj(r, z, phi)
    return conj(central_fdm(3,1)(z1 -> Ez(r, z1, phi), z))
end

function Fz(r, z, phi)
    term1 = 2 * pi / eta_var
    term2 = alpha_var * (Ex(r,z) * dEx_conj(r,z) + Ez(r,z,phi) * dEz_conj(r,z,phi))
    return term1 * real(term2)
end

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Fzvalues = []
global count = 0
for z in zvalues
    push!(Fzvalues, Fz(0, -z, 0))
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

save("Dipole_force.png", fig)
