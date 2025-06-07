using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, Interpolations

max = 3
tam = 50
zs = LinRange(-max, max, tam)
w0 = 0.3
z0 = 1
intensity1 = [w0 * sqrt(1 + (z/z0)^2) for z in zs]
intensity2 = [-w0 * sqrt(1 + (z/z0)^2) for z in zs]

fig = Figure(fontsize=30, px_per_unit = 2, size= 3 .*(420, 170), figure_padding = 25)
ax = Axis(fig[1, 1], xlabel = L"z", 
    ylabel = L"W(z)", aspect=DataAspect(), yticksvisible=false, yticklabelsvisible=false,
    xticks = (-3:1:3, [L"-3z_0", L"-2z_0", L"-1z_0", L"0", L"1z_0", L"2z_0", L"3z_0"]))

#lines!(ax, zs, intensity1, color="black")
#lines!(ax, zs, intensity2, color="black")

points = [Point2f[(zs[i], intensity1[i]) for i in range(1, tam)]; Point2f[(zs[tam + 1 - i], intensity2[tam + 1 - i]) for i in range(1, tam)]]
poly!(ax, points, color="red", alpha=0.2)

zcircs = LinRange(0.25,max,12)
for z in zcircs
    radius = z * (1 + (z0 / z)^2)
    theta = 0
    while radius * sin(theta) < w0 * sqrt(1 + ((z - radius * (1 - cos(theta))) / z0)^2)
        theta += 0.005
    end
    angles = LinRange(-theta, theta, tam)
    wz = w0 * sqrt(1 + (z / z0)^2)

    local rcirc = radius .* sin.(angles)
    local zcirc = z .- radius .* (1 .- cos.(angles))
    lines!(ax, zcirc, rcirc, color="red")
    lines!(ax, -zcirc, rcirc, color="red")
end


save("raio_feixe_gauss.png", fig)