using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions

tam = 500
xs = LinRange(-8, 8, tam)
ys = LinRange(-8, 8, tam)
zs = LinRange(-50, 50, tam)
kr = 0.8
intensity_bessel = [(besselj0(kr * x)).^2 for z in zs, x in xs]

fig = Figure(fontsize=36, figure_padding=20)
ax1 = Axis(fig[1, 1], xlabel = L"z", ylabel = L"y",
    aspect=2.5, xticklabelsvisible=false, yticklabelsvisible=false)
hm1 = heatmap!(ax1, zs, xs, abs.(intensity_bessel))

intensity_bessel2 = [(besselj0(kr * sqrt(x^2 + y^2))).^2 for x in xs, y in ys]

ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y",
    aspect=1, xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
    ylabelvisible=false)
hm2 = heatmap!(ax2, xs, ys, abs.(intensity_bessel2))

colgap!(fig.layout, Relative(0.02))
colsize!(fig.layout, 1, Aspect(1, 2.5))
colsize!(fig.layout, 2, Aspect(1, 1))
#rowgap!(fig.layout, 0.01)
resize_to_layout!(fig)

save("feixe_bessel.png", fig)
Makie.inline!(true)
current_figure()
