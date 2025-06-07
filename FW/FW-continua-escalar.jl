using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions

omega = 3.14 * 10^(15)
c = 299792458
K = 2 * omega / c

function F1(x, Z)
    if(x < 0)
        return 0
    end
    a = 3/Z
    return exp(a*(x - Z))
end

function F2(x, Z)
    q = 2
    return exp(-q*(x / Z)^2)
end

function F3(x, Z)
    q = 2
    return exp(-q*(x / Z)^8)
end

function F4(x, Z)
    q = 1.6 * 10^6
    return besselj0(q * x)
end

function FW(func, N, r, z, Z)
    terms = range(-N, N)
    sum_mackinnon = 0 + 0im
    for n in terms
        sum_mackinnon += func(- 2 * n * pi / K, Z) * sinc(1/pi * sqrt((omega * r/c)^2 + (z * omega / c + pi * n)^2))
    end
    return sum_mackinnon
end

tam = 500
R = 2 * 10^(-6)

Z1 = 10 * 10^(-6)
N1 = Int(round(Z1 * omega / (pi * c), RoundDown))
r1 = LinRange(-R, R, tam)
z1 = LinRange(-Z1, Z1, tam)
intensity1 = [abs2(FW(F1, N1, r, z, Z1)) for r in r1, z in z1]

Z2 = 1.6 * 10^(-6)
N2 = Int(round(Z2 * omega / (pi * c), RoundDown))
r2 = LinRange(-R, R, tam)
z2 = LinRange(-Z2, Z2, tam)
intensity2 = [abs2(FW(F2, N2, r, z, Z2)) for r in r2, z in z2]

Z3 = 2 * 10^(-6)
N3 = Int(round(Z3 * omega / (pi * c), RoundDown))
r3 = LinRange(-R, R, tam)
z3 = LinRange(-Z3, Z3, tam)
intensity3 = [abs2(FW(F3, N3, r, z, Z3)) for r in r3, z in z3]

Z4 = 15 * 10^(-6)
N4 = Int(round(Z4 * omega / (pi * c), RoundDown))
r4 = LinRange(-R, R, tam)
z4 = LinRange(-Z4, Z4, tam)
intensity4 = [abs2(FW(F4, N4, r, z, Z4)) for r in r4, z in z4]

# max_intensity = maximum(abs.(intensity))
# intensity = intensity / max_intensity

fig = Figure(figure_padding=18)
ax1 = Axis3(fig[1, 1], 
    azimuth = -0.25 * pi, elevation = 0.12 * pi, xlabel = L"\rho \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"|\Psi|^2")
ax2 = Axis3(fig[1, 2], 
    azimuth = -0.25 * pi, elevation = 0.12 * pi, xlabel = L"\rho \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"|\Psi|^2")
ax3 = Axis3(fig[2, 1], 
    azimuth = -0.25 * pi, elevation = 0.12 * pi, xlabel = L"\rho \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"|\Psi|^2")
ax4 = Axis3(fig[2, 2], 
    azimuth = -0.25 * pi, elevation = 0.12 * pi, xlabel = L"\rho \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"|\Psi|^2")

Label(fig[1,1][1, 1, TopLeft()], "(a)", fontsize = 16,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[1,2][1, 1, TopLeft()], "(b)", fontsize = 16,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[2,1][1, 1, TopLeft()], "(c)", fontsize = 16,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[2,2][1, 1, TopLeft()], "(d)", fontsize = 16,
    padding = (0, 10, -5, 0), halign = :right)

r1 *= 10^6
r2 *= 10^6
r3 *= 10^6
r4 *= 10^6
z1 *= 10^6
z2 *= 10^6
z3 *= 10^6
z4 *= 10^6

surface!(ax1, r1, z1, intensity1)
surface!(ax2, r2, z2, intensity2)
surface!(ax3, r3, z3, intensity3)
surface!(ax4, r4, z4, intensity4)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("FW_continua_escalar.png", fig)
