using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes

tam = 1000
xs = LinRange(-10, 10, tam)
zs = LinRange(5, 50, tam)
k = 7
z0 = 1im * 10
intensity = [1/z * exp(-1im * k * x^2 / (2 * z)) for x in xs, z in zs]
intensity2 = [1/(z + z0) * exp(-1im * k * x^2 / (2 * (z + z0))) for x in xs, z in zs]

fig = Figure()
ax = Axis3(fig[1, 1], title="Distribuição de intensidade do feixe",
    azimuth = 0.45 * pi, xlabel = "x", ylabel = "z", zlabel = "Intensidade")
ax2 = Axis3(fig[1, 2], title="Distribuição de intensidade do feixe",
    azimuth = 0.45 * pi, xlabel = "x", ylabel = "z", zlabel = "Intensidade")
surface!(ax, xs, zs, real.(intensity))
surface!(ax2, xs, zs, abs.(intensity2))

save("feixe_gauss.png", fig)
