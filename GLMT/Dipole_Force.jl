using GLMakie, ColorSchemes, ColorTypes, Base.Threads
using MakieThemes, SpecialFunctions, Interpolations, Integrals, FiniteDifferences, Arblib

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
Q = 0.85 * k
E0 = 1
eta_var = sqrt(mu / epsilon)
mult_factor = BigFloat("1e+24")

nmax = Int(round((Zmax + z0) * k / pi, RoundUp))
nmin = Int(round((-Zmax + z0) * k / pi, RoundDown))

tam = 500

# calculo para a << lambda
#np = BigFloat(1.2)
#size = k * BigFloat(10^-9)
# errado a partir daqui, faltando algo
#alpha_var = (np^2 - 1) / (np^2 + 2) * size^3 / k^3 - 1im * 2 / 3 * ((np^2 - 1) / (np^2 + 2))^2 * size^6 / k^3


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
end

function ricatti_besselh(n, x)
    return x * custom_hankel2(n, x)
end

function dricatti_besselh(n, x)
    return forward_fdm(5,1)(x1 -> ricatti_besselh(n, x1), x)
    #return - n * custom_hankel2(n, x) + x * custom_hankel2(n-1, x)
end

function an(n)
    term1 = musp * ricatti_besselj(n, alpha) * dricatti_besselj(n, betah) - 
        mu * M * ricatti_besselj(n, betah) * dricatti_besselj(n, alpha)
    term2 = musp * ricatti_besselh(n, alpha) * dricatti_besselj(n, betah) - 
        mu * M * ricatti_besselj(n, betah) * dricatti_besselh(n, alpha)
    return term1 / term2
end
a1_mie = an(1)
alpha_var = - 1im * 3 / (2 * k^3) * a1_mie
println(a1_mie)


@inline function F(z)
    term1 = exp(1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

function g(z)
    prefactor = (1/4) * sqrt(big(pi) / 2) * Zmax
    arg_erf = Acb((2 * sqrt(2) * z) / Zmax + (im * Q * Zmax) / (4 * sqrt(2)))

    return prefactor * exp(- (Q^2 * Zmax^2) / 32) * Complex{BigFloat}(erf(arg_erf))
end

function Ex(r, z)
    values = 0
    for n in range(nmin,nmax)
        zf = - big(pi) * n / k
        values += F(zf) * sinc.(1/big(pi) * sqrt((k * r)^2 + (z * k + big(pi) * n)^2))
    end
    return E0 * values
end

function Ez(r, z, phi)
    values = 0
    for n in range(nmin,nmax)
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
    for n in range(nmin,nmax)
        zf = - big(pi) * n / k
        values += F(zf) * dsinc_analytic(r, z, n)
    end
    return conj(E0 * values)
end

function dEz_conj(r, z, phi)
    return conj(central_fdm(3,1)(z1 -> Ez(r, z1, phi), z))
end

function Fz(r, z, phi)
    term1 = 2 * big(pi) / eta_var
    term2 = alpha_var * (Ex(r,z) * dEx_conj(r,z) + Ez(r,z,phi) * dEz_conj(r,z,phi))
    return term1 * real(term2) * mult_factor
end

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")

Fzvalues = Array{BigFloat}(undef, tam)
Threads.@threads for count in 1:tam
    global Fzvalues[count] =  Fz(0, -zvalues[count], 0)
    if count % 100 == 0
        println(count, " pontos calculados")
    end
end

# plots
fig = Figure(figure_padding=40, size= 3 .*(500, 200), fontsize = 25)
ax1 = Axis(fig[1, 1], title = "(a)",
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.5)
ax2 = Axis(fig[1, 2], title = "(b)", yticklabelsvisible=false, yticksvisible=false,
    xlabel = L"z_0 \, (\mu m)", ylabel = L"F_z \, (pm^2)", aspect=1.5)
ax_inset = Axis(fig[1, 1],
    width=Relative(0.2),
    height=Relative(0.2),
    halign=0.05, valign=0.85,
    xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
    ylabelvisible=false, xticksvisible=false)

function FW(r, z)
    sum_mackinnon = 0 + 0im
    for n in range(nmin, nmax)
        sum_mackinnon += F(- n * pi / k) * sinc(1/pi * sqrt((k * r)^2 + (z * k + pi * n)^2))
    end
    return sum_mackinnon
end

Fterms = F.(- big(pi) .* range(nmin,nmax) ./ k)
function S(x)
    return 1/K * reduce(+, Fterms .* exp.(1im * 2 * pi * range(nmin,nmax) * x / K))
end

Fplot = abs2.(F.(zvalues))
intensityFW = [abs2(FW(0, z)) for z in zvalues]
spectra = [abs2(S(x)) for x in range(-k,k,tam)]

zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="black", label=L"|F(z)|^2")
lines!(ax1, zvalues, intensityFW, color="blue", label=L"|\Psi(0,z)|^2", linestyle=:dash)
lines!(ax_inset, zvalues, spectra, color="black")
lines!(ax2, zvalues, Fzvalues, color="black")

axislegend(ax1, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("Dipole_force.png", fig)
