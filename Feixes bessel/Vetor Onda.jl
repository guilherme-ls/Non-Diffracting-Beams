using GLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations, LinearAlgebra

fig = Figure(size= 2 .*(500, 500), fontsize=24, figure_padding=2)
ax = Axis3(fig[1, 1], azimuth = 0.25 * pi, elevation = 0.12 * pi,
    xlabel = L"x", ylabel = L"y", zlabel = L"z", aspect=(1, 1, 1))

hidedecorations!(ax)
hidespines!(ax)

angle_finish = LinRange(0, 2*pi, 10)
origin = [Point3f([0, 0, 0]) for a in angle_finish]
finish = [Vec3f([cos(a), sin(a), 1]) for a in angle_finish]
angle_circle = LinRange(0, 2*pi, 100)
circling = [Point3f([cos(a), sin(a), 1]) for a in angle_circle]

size = Vec3f(0.045, 0.045, 0.08)
lines!(ax, circling, color=:black, linewidth=3)
arrows!(ax, origin, finish .* (norm.(finish).-size[3])./norm.(finish),
    arrowsize=size, linewidth=0.01)

origin_coord = [Point3f([0,0,0]) for i in range(1,3)]
end_coord= [Vec3f([0,0,2]), Vec3f([1.2,0,0]), Vec3f([0,1.2,0])]
arrows!(ax, origin_coord, end_coord .* (norm.(end_coord).-size[3])./norm.(end_coord),
    arrowsize=size, linewidth=0.015, color=:gray)

a_circle = LinRange(0, pi/4, 100)
dist = 1.8
a_circle_pos = [Point3f([-sin(a) * sqrt(2)/2 * dist, sin(a) * sqrt(2)/2 * dist, cos(a) * dist]) for a in a_circle]
angle_line = [Point3f([0,0,0]), Point3f(1.1 .*[-dist/2,dist/2,dist/sqrt(2)])]

lines!(ax, a_circle_pos, color=:gray, linewidth=2, linestyle=:dash)
lines!(ax, angle_line, color=:gray, linewidth=2)

text!(ax, -0.02, 0.02, 1.8, text=L"z")
text!(ax, 1.1, -0.04, 0, text=L"x")
text!(ax, -0.04, 1.05, 0, text=L"y")
text!(ax, -0.5, 0.5, 1.7, text=L"\theta")

resize_to_layout!(fig)

save("vetor_onda_bessel.png", fig)
Makie.inline!(true)
current_figure()