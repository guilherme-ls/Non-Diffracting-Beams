using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions

tam = 500
xs = LinRange(-8, 8, tam)
ys = LinRange(-8, 8, tam)
zs = LinRange(-50, 50, tam)
k = 10
kr = 0.8
z0 = 1im * 20
intensity_gauss = [(1 * exp(-1im * k * x^2 / (2 * (z + z0)))).^2 for z in zs, x in xs]
intensity_bessel = [(besselj0(kr * x)).^2 for z in zs, x in xs]

fig = Figure(fontsize=18, figure_padding=10)
ax = Axis(fig[1, 1], xlabel = L"z", ylabel = L"y",
    aspect=2.5, xticklabelsvisible=false, yticklabelsvisible=false)
ax2 = Axis(fig[2, 1], xlabel = L"z", ylabel = L"y",
    aspect=2.5, xticklabelsvisible=false, yticklabelsvisible=false)
hm1 = heatmap!(ax, zs, xs, abs.(intensity_gauss))
hm2 = heatmap!(ax2, zs, xs, abs.(intensity_bessel))

intensity_gauss2 = [(1 * exp(-1im * k * (x^2 + y^2) / (2 * (30 + z0)))).^2 for x in xs, y in ys]
intensity_bessel2 = [(besselj0(kr * sqrt(x^2 + y^2))).^2 for x in xs, y in ys]

ax3 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y",
    aspect=1, xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
    ylabelvisible=false)
ax4 = Axis(fig[2, 2], xlabel = L"x", ylabel = L"y",
    aspect=1, xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
    ylabelvisible=false)
hm1 = heatmap!(ax3, xs, ys, abs.(intensity_gauss2))
hm2 = heatmap!(ax4, xs, ys, abs.(intensity_bessel2))

Label(fig[1,1][1, 1, TopLeft()], "a)", fontsize = 22, font = :bold,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[2,1][1, 1, TopLeft()], "b)", fontsize = 22, font = :bold,
    padding = (0, 10, -5, 0), halign = :right)

colgap!(fig.layout, Relative(0.02))
colsize!(fig.layout, 1, Aspect(1, 2.5))
colsize!(fig.layout, 2, Aspect(1, 1))
#rowgap!(fig.layout, 0.01)
resize_to_layout!(fig)

save("feixe_bessel_gauss.png", fig)
Makie.inline!(true)
current_figure()
