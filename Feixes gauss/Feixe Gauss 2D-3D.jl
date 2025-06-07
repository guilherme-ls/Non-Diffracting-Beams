using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations

tam = 500
xs = LinRange(-8, 8, tam)
zs = LinRange(-50, 50, tam)
k = 15
z0 = 1im * 12
intensity = [(1/(z + z0) * exp(-1im * k * x^2 / (2 * (z + z0)))).^2 for x in xs, z in zs]
max_intensity = maximum(abs.(intensity))
intensity = intensity / max_intensity

fig = Figure(figure_padding=2)
ax = Axis(fig[1, 1],
    xlabel = L"x", ylabel = L"z", aspect=0.4)
ax2 = Axis3(fig[1, 2:3], 
    azimuth = 0.35 * pi, elevation = 0.12 * pi, xlabel = L"x", ylabel = L"z", 
    zlabel = L"I/I_0", aspect=(1, 3, 1))

hm = heatmap!(ax, xs, zs, abs.(intensity))
surface!(ax2, xs, zs, abs.(intensity))
contour3d!(ax2, xs, zs, abs.(intensity), 
    levels=(minimum(abs.(intensity)):0.1:1), linewidth=1.5)

colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 0.01)
#resize_to_layout!(fig)

save("feixe_gauss_2D_3D.png", fig)
