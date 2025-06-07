using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations

max = 3
tam = 100
zs = LinRange(-max, max, tam)
w0 = 1
z0 = 1
intensity1 = [w0 * sqrt(1 + (z/z0)^2) for z in zs]
intensity2 = [w0 * (-1) * sqrt(1 + (z/z0)^2) for z in zs]

assintotax = [-max, max]
assintotay = [-max * w0 / z0, max * w0 / z0]

fig = Figure(fontsize=36, px_per_unit = 2, size= 3 .*(500, 170), figure_padding = 25)
ax = Axis(fig[1, 1], xlabel = L"z", 
    ylabel = L"W(z)", aspect=4, limits=(-max, max, 0, max+0.2),
    xticks = (-max:1:max, [L"-3z_0", L"-2z_0", L"-1z_0", L"0", L"1z_0", L"2z_0", L"3z_0"]),
    yticks = (0:1:max, [L"0", L"\omega_0", L"2 \omega_0", L"3 \omega_0"]))

lines!(ax, zs, intensity1, color="black")
#lines!(ax, zs, intensity2, color="black")
lines!(ax, assintotax, assintotay, color="gray", linestyle=:dash)
lines!(ax, assintotax, -assintotay, color="gray", linestyle=:dash)

#colsize!(fig.layout, 1, Aspect(1, 4.0))

save("limites_feixe_gauss.png", fig)