using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions

tam = 10000
zmax = 5
theta = 0.01
radius = zmax * tan(theta)
xs = LinRange(-radius, radius, tam)
angle = LinRange(0, 2 * pi, tam)
lambda = 430 * 10^(-9)
k = 2 * pi / lambda
kr = k * sin(theta)
A = 1 / tam

fig = Figure(fontsize=18, figure_padding=30, size = (800, 800))

plane = [0.0 + 0.0im for x in xs]
xi = 1
for x in xs
    global xi
    for a in angle
        plane[xi] += A * exp(1im * kr * x * cos(a))
    end
    xi += 1
end

ax = Axis(fig[1, 1], xlabel = L"\rho", ylabel = L"I", aspect=1)
lines!(ax, xs, abs.(plane))

#rowsize!(fig.layout, 1, Aspect(1, 2))
#rowgap!(fig.layout, 0.01)
#resize_to_layout!(fig)

save("feixe_bessel_axicon.png", fig)
Makie.inline!(false)
current_figure()