using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations, SpecialFunctions

# constantes
L = 0.5
l1 = L/10
l2 = 3*L/10
l3 = 4*L/10
l4 = 6*L/10
l5 = 7*L/10
l6 = 9*L/10

c = 299792458
w = 2.98 * 10^15
Q = 0.9998 * w / c
N = 20

# função a ser aproximada
function funcF(z)
    if l1 <= z && z <= l2
        return -4*(z - l1)*(z - l2)/((l2 - l1)^2)
    end
    if l3 <= z && z <= l4
        return 1
    end
    if l5 <= z && z <= l6
        return -4*(z - l5)*(z - l6)/((l6 - l5)^2)
    end
    return 0
end

N = 25
tam = 500
zvalues = LinRange(0, L, tam)
nvalues = range(-N, N)

# funções para os parâmetros A e B
function An(n)
    return 1/tam * reduce(+, [funcF(z) * exp(-1im * 2 * pi * n * z /L) for z in zvalues])
end
Avalues = [An(n) for n in nvalues]

function Bn(n)
    return Q + 2 * pi * n / L
end

function kr(n)
    return sqrt(w^2 - Bn(n)^2)
end
Kvalues = [kr(n) for n in nvalues]

function FW(r, z, t)
    sum_bessel = reduce(+, Avalues .* besselj0.(Kvalues * r) .* exp.(1im * 2 * pi * nvalues * z / L))
    return exp(- 1im * w * t + 1im * Q * z) * sum_bessel
end

r1 = LinRange(-5*10^-4, 5*10^-4, tam + 1)
intensity = [abs2(FW(r, z, 0)) for r in r1, z in zvalues]
intensity_lines = [abs2(FW(0, z, 0)) for z in zvalues]

# max_intensity = maximum(abs.(intensity))
# intensity = intensity / max_intensity

r1 *= 1000

fig = Figure(figure_padding=50, size= 2 .*(500, 250), fontsize = 16)
ax = Axis(fig[1, 1],
    xlabel = L"z", ylabel = L"|\Psi|^2", aspect=1.3)
ax2 = Axis3(fig[1, 2], 
    azimuth = -0.25 * pi, elevation = 0.12 * pi, xlabel = L"\rho \, (mm)",
    ylabel = L"z \, (m)", zlabel = L"|\Psi|^2")

lines!(ax, zvalues, intensity_lines, label=L"|\Psi_{FW}|^2")
lines!(ax, zvalues, abs2.(funcF.(zvalues)), color="black", label=L"|F(z)|^2")
surface!(ax2, r1, zvalues, intensity)

Label(fig[1,1][1, 1, TopLeft()], "(a)", fontsize = 18,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[1,2][1, 1, TopLeft()], "(b)", fontsize = 18,
    padding = (0, 10, -5, 0), halign = :right)

axislegend(ax, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowsize!(fig.layout, 1, Relative(0.6))
#trim!(fig.layout)
resize_to_layout!(fig)

save("FW_discreta.png", fig)
