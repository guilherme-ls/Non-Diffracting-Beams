using WGLMakie, ColorSchemes, ColorTypes
using MakieThemes, SpecialFunctions, Interpolations, Integrals

lambda = 532 * 10^-9
c = 299792458
omega = 2 * pi * c / lambda
K = 2 * omega / c
Z = 12.61 * 10^-6
#N = Int(round(Z * omega / (pi * c), RoundDown))
N = 52
Q = 0.7 * omega / c
E0 = 1

tam = 500

function F(z)
    term1 = exp(1im * Q * z)
    term2 = exp(-0.5 * (z/Z)^8)
    term3 = airyai(3 * pi * (z/Z - 0.8))
    return term1 * term2 * term3
end

function f(z, p)
    return F(z)
end

terms = range(-N, N)
Fvalues = F.(- 2 * terms * pi / K)

function g(n)
    prob = IntegralProblem(f, (-1.1 * Z, -2 * pi * n /K))
    return solve(prob, HCubatureJL(), reltol = 1e-20, abstol = 1e-20).u
end

function gsum(n)
    values = LinRange(-10 * Z, -2 * pi * n /K, 100000)
    soma = reduce(+, F.(values) * (values[2] - values[1]))
    return soma
end
gvalues = gsum.(terms)
#gvalues = g.(terms)
#gvalues = F.(- 2 * terms * pi / K)

function Ex(r, z)
    values = Fvalues .* sinc.(1/pi .* sqrt.((omega .* r./c).^2 .+ (z .* omega ./ c .+ pi .* terms).^2))
    sum_mackinnon = reduce(+, values)
    return E0 * sum_mackinnon
end

function Ez(r, z, phi)
    if r == 0
        return 0
    end
    values = gvalues .* omega^2 .* r ./ c^2 .* (cos.(sqrt.((omega .* r./c).^2 .+ (z .* omega ./ c .+ pi .* terms).^2))./((omega * r/c)^2 .+ (z * omega / c .+ pi .* terms).^2) - sin.(sqrt.((omega * r/c)^2 .+ (z * omega / c .+ pi .* terms).^2))./((omega * r/c)^2 .+ (z * omega / c .+ pi .* terms).^2).^1.5)
    sum_mackinnon = reduce(+, values)
    return -E0 * cos(phi) * sum_mackinnon
end

function S(k)
    return 1/K * reduce(+, Fvalues .* exp.(1im * 2 * pi * terms * k / K))
end

R = 1 * 10^(-6)

rvalues = LinRange(-R, R, tam)
zvalues = LinRange(-Z, Z, tam)
kvalues = LinRange(-omega/c, omega/c, tam)
intensity = [abs2(Ex(0, z)./E0) for z in zvalues]
spectra = [abs(S(k)) for k in kvalues]

# unit vectors:
xuv, yuv, zuv = Point3f(1,0,0), Point3f(0,1,0), Point3f(0,0,1)

# transform x,y,z coordinates to points in 3D:
points_plane1 = [ix * xuv + iz * yuv + 0 * zuv for ix in rvalues for iz in zvalues]
points_plane2 = [0 * xuv + iz * yuv + iy * zuv for iy in rvalues for iz in zvalues]
points_plane3 = [ix * xuv + 9 * 10^-6 * yuv + iy * zuv for ix in rvalues for iy in rvalues]

# find values of F in all in-plane points:
fcolors11 = [abs2(Ex(x, z)) for x in rvalues for z in zvalues]
fcolors21 = [abs2(Ex(y, z)) for y in rvalues for z in zvalues]
fcolors31 = [abs2(Ex(sqrt(x^2 + y^2), 9 * 10^-6)) for x in rvalues for y in rvalues]

# find values of F in all in-plane points:
fcolors12 = [abs2(Ez(x, z, 0)) for x in rvalues for z in zvalues]
fcolors22 = [abs2(Ez(y, z, pi/2)) for y in rvalues for z in zvalues]
fcolors32 = [abs2(Ez(sqrt(x^2 + y^2), 9 * 10^-6, atan(y,x))) for x in rvalues for y in rvalues]

# plots
fig = Figure(figure_padding=40, size= 3 .*(500, 400), fontsize = 20)
ax1 = Axis(fig[1, 1], 
    xlabel = L"z \, (\mu m)", ylabel = L"|\Psi|^2")
ax2 = Axis(fig[1, 2], 
    xlabel = L"k_z / k", ylabel = L"|S(k_z)|")
ax3 = Axis3(fig[2, 1], 
    azimuth = -0.15 * pi, elevation = 0.14 * pi, xlabel = L"x \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"y \, (\mu m)")
ax4 = Axis3(fig[2, 2], 
    azimuth = -0.15 * pi, elevation = 0.14 * pi, xlabel = L"x \, (\mu m)",
    ylabel = L"z \, (\mu m)", zlabel = L"y \, (\mu m)")

Label(fig[1,1][1, 1, TopLeft()], "(a)", fontsize = 22,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[1,2][1, 1, TopLeft()], "(b)", fontsize = 22,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[2,1][1, 1, TopLeft()], "(c)", fontsize = 22,
    padding = (0, 10, -5, 0), halign = :right)
Label(fig[2,2][1, 1, TopLeft()], "(d)", fontsize = 22,
    padding = (0, 10, -5, 0), halign = :right)

Fplot = abs2.(F.(zvalues))

rvalues *= 10^6
zvalues *= 10^6
kvalues *= c / omega

lines!(ax1, zvalues, intensity, label=L"|\Psi_{FWC}|^2", color="blue")
lines!(ax1, zvalues, Fplot, color="black", label=L"|F(z)|^2")
lines!(ax2, kvalues, spectra, color="black")

append!(points_plane1, points_plane2, points_plane3)
append!(fcolors11, fcolors21, fcolors31)
append!(fcolors12, fcolors22, fcolors32)

points_plane1 *= 10^6
points_plane2 *= 10^6
points_plane3 *= 10^6

scatter!(ax3, points_plane1; alpha=0.7,
    color=fcolors11, colormap=:viridis, markersize=0.1)

scatter!(ax4, points_plane1; alpha=0.7,
    color=fcolors12, colormap=:viridis, markersize=0.1)

println(maximum(fcolors12))
println(minimum(fcolors12))
println(maximum(abs2.(gvalues)))
println(minimum(abs2.(gvalues)))

axislegend(ax1, position=:lt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("FW_continua_vetorial.png", fig)
