using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations

tam = 1000
xs = LinRange(-10, 10, tam)
zs = LinRange(5, 50, tam)
zs2 = LinRange(-50, 50, tam)
k = 7
z0 = 1im * 10
intensity = [1/z * exp(-1im * k * x^2 / (2 * z)) for x in xs, z in zs]
intensity2 = [1/(z + z0) * exp(-1im * k * x^2 / (2 * (z + z0))) for x in xs, z in zs2]

fig = Figure()
ax = Axis(fig[1, 1], title="Feixe parabolóide",
    xlabel = "x", ylabel = "z")
ax2 = Axis(fig[1, 2], title="Feixe gaussiano",
    xlabel = "x", ylabel = "z")

heatmap!(ax, xs, zs, real.(intensity))
heatmap!(ax2, xs, zs, abs.(intensity2))

save("feixe_gauss_2D.png", fig)
