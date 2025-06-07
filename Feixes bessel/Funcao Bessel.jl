using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions

max = 20
tam = 100
xs = LinRange(0, max, tam)
xs2 = LinRange(0.1, max, tam)

fig = Figure(fontsize=18)
ax = Axis(fig[1, 1], xlabel = L"x", ylabel = L"y", aspect=1.4, limits=(0, max, -0.5, 1.1))

bessel10 = besselj.(0, xs)
bessel11 = besselj.(1, xs)
bessel12 = besselj.(2, xs)
lines!(ax, xs, bessel10, label=L"J_0(x)", color="red")
lines!(ax, xs, bessel11, label=L"J_1(x)", color="green")
lines!(ax, xs, bessel12, label=L"J_2(x)", color="blue")

bx = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y", aspect=1.4, limits=(0, max, -1.5, 1))
bessel20 = bessely.(0, xs2)
bessel21 = bessely.(1, xs2)
bessel22 = bessely.(2, xs2)
lines!(bx, xs2, bessel20, label=L"Y_0(x)", color="red")
lines!(bx, xs2, bessel21, label=L"Y_1(x)", color="green")
lines!(bx, xs2, bessel22, label=L"Y_2(x)", color="blue")

colsize!(fig.layout, 1, Aspect(1, 1.4))
colsize!(fig.layout, 2, Aspect(1, 1.4))
resize_to_layout!(fig)
axislegend(ax)
axislegend(bx, position=:rb)

save("funcao_bessel.png", fig)
Makie.inline!(true)
current_figure()