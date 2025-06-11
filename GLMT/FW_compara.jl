using WGLMakie, ColorSchemes, ColorTypes, JLD2
using MakieThemes, SpecialFunctions, HypergeometricFunctions, AssociatedLegendrePolynomials

lambda = 1064 * 10^-9
c = 299792458
omega = 2 * pi * c / lambda
k = omega / c
K = 2 * k
Zmax = 22 * 10^-6
pmax = 1000
nmax = 47
z0 = 0
Q = 0.9 * k
E0 = 1

lmax = Int(round(Zmax * k / pi, RoundUp))
lmin = Int(round(-Zmax * k / pi, RoundDown))

tam = 100

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

function F(z)
    term1 = exp(- 1im * Q * z)
    term2 = exp(- 8 * (z/Zmax)^2)
    return term1 * term2
end

Fterms = F.(pi .* range(-lmax,lmax) ./ k)

function gTM(n, m)
    term1 = - (1.0im)^m / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(n-m+1) - loggamma(n + abs(m) + 1)) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint(pi * (l + p) + k * z0)
        end
        term2 += (B1_data[(n, p)] * (m + 1) + B2_data[(n, p)] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function gTE(n, m)
    term1 = (1.0im)^(m + 1) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(n-m+1) - loggamma(n + abs(m) + 1)) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pmax,pmax)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint(pi * (l + p) + k * z0)
        end
        term2 += (B3_data[(n, p)] * (m + 1) - B4_data[(n, p)] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function cn(n)
    return 1 / (1.0im * k) * (-1.0im)^n * (2*n+1) / (n * (n+1))
end

cterms = cn.(range(1,nmax))

# calculates Legendre associated function and derivatives
dP, P = dlegendre(Val(:unnormalized), 1, nmax, 1; ph_term=true)

#println(dP)

function ErTM(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 +=  gTM(n,m) * Plm(n, 1, 1) * exp(1im * m * phi)
        end
        term1 += cterms[n] * n * (n+1) / r * sphericalbesselj(n, k * r) * term2
    end
    return E0 * term1
end

function EthetaTM(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 +=  gTM(n,m) * (Plm(n, 1, cos(theta + 0.001)) - Plm(n, 1, cos(theta)))/0.001 * exp(1im * m * phi)
        end
        term1 += cterms[n] * (k * (r + 0.001) * sphericalbesselj(n, k * (r + 0.001)) - k * r * sphericalbesselj(n, k * r))/ 0.001 * term2
    end
    return E0 * term1 / r
end

function EphiTM(r, theta, phi)
    if theta == 0
        return 0
    end
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += m * gTM(n,m) * P[n + 1, 2] / sin(theta) * exp(1im * m * phi)
        end
        term1 += cterms[n] * (k * sphericalbesselj(n-1, k * r) - n / r * sphericalbesselj(n, k * r)) * term2
    end
    return E0 * term1 * 1im
end

function ErTE(r, theta, phi)
    return 0
end

function EthetaTE(r, theta, phi)
    if theta == 0
        return 0
    end
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += m * gTE(n,m) * P[n + 1, 2] / sin(theta) * exp(1im * m * phi)
        end
        term1 += cterms[n] * sphericalbesselj(n, k * r) * term2
    end
    return E0 * term1 * k
end

function EphiTE(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += gTE(n,m) * dP[n + 1, 2] * exp(1im * m * phi)
        end
        term1 += cterms[n] * sphericalbesselj(n, k * r) * term2
    end
    return E0 * term1 * 1im * k
end

function Er(r, theta, phi)
    return ErTE(r, theta, phi) + ErTM(r, theta, phi)
end

function Etheta(r, theta, phi)
    return EthetaTE(r, theta, phi) + EthetaTM(r, theta, phi)
end

function Ephi(r, theta, phi)
    return EphiTE(r, theta, phi) + EphiTM(r, theta, phi)
end

function Ex(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    theta = 0 #acos(z/r)
    phi = 0 #atan(y/x)
    comp1 = 0 # Er(r, theta, phi) * sin(theta) * cos(phi)
    comp2 = Etheta(r, theta, phi) * cos(theta) * cos(phi)
    comp3 = 0 #- Ephi(r, theta, phi) * sin(phi)
    return comp1 + comp2 + comp3
end

function Ey(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    theta = 0 #acos(z/r)
    phi = 0 #atan(y/x)
    comp1 = 0#Er(r, theta, phi) * sin(theta) * sin(phi)
    comp2 = 0#Etheta(r, theta, phi) * cos(theta) * sin(phi)
    comp3 = Ephi(r, theta, phi) * cos(phi)
    return comp1 + comp2 + comp3
end

function Ez(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    theta = acos(z/r)
    phi = 0 #atan(y/x)
    comp1 = Er(r, theta, phi) * cos(theta)
    comp2 = 0 # - Etheta(r, theta, phi) * sin(theta)
    return comp1 + comp2
end

function FW(r, z)
    sum_mackinnon = 0 + 0im
    for n in range(-lmax, lmax)
        sum_mackinnon += F(- 2 * n * pi / K) * sinc(1/pi * sqrt((omega * r/c)^2 + (z * omega / c + pi * n)^2))
    end
    return sum_mackinnon
end

R = 1 * 10^(-6)

zvalues = LinRange(-Zmax, Zmax, tam)
println("Starting Calculation")
intensityX = [abs2(Ex(0, 0, z)) for z in zvalues]
intensityZ = [abs2(Ez(0, 0, z)) for z in zvalues]
intensityFW = [abs2(FW(0, z)) for z in zvalues]

# normaliza, remover depois
val_normalize = maximum(intensityX)
intensityX = intensityX ./ val_normalize
#intensityZ = intensityZ ./ maximum(intensityZ)

# plots
fig = Figure(figure_padding=40, size= 3 .*(600, 200), fontsize = 25)
ax1 = Axis(fig[1, 1], 
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.3)
ax2 = Axis(fig[1, 2], 
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.3)
ax3 = Axis(fig[1, 3], 
    xlabel = L"z \, (\mu m)", ylabel = L"ln(erro)", aspect=1.3)

Label(fig[1,1][1, 1, TopLeft()], "(a)", fontsize = 25,
    padding = (0, 10, -7, 0), halign = :right)
Label(fig[1,2][1, 1, TopLeft()], "(b)", fontsize = 25,
    padding = (0, 10, -7, 0), halign = :right)
Label(fig[1,3][1, 1, TopLeft()], "(c)", fontsize = 25,
    padding = (0, 10, -7, 0), halign = :right)

Fplot = abs2.(F.(zvalues))

zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="black", label=L"|F(z)|^2")
lines!(ax1, zvalues, intensityFW, color="blue", label=L"|\Psi(0,z)|^2", linestyle=:dash)
lines!(ax2, zvalues, intensityFW, color="black", label=L"|\Psi(0,z)|^2")
lines!(ax2, zvalues, intensityX, label=L"|E_x|^2", color="blue", linestyle=:dash)
lines!(ax2, zvalues, intensityZ, label=L"|E_z|^2", color="red", linestyle=:dashdot)
lines!(ax3, zvalues, log.(abs.(intensityX .- intensityFW)./intensityFW), color="black")

axislegend(ax1, position=:rt)
axislegend(ax2, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("FW_compara.png", fig)
