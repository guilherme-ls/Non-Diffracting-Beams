using GLMakie, ColorSchemes, ColorTypes, JLD2, FiniteDifferences
using MakieThemes, SpecialFunctions, AssociatedLegendrePolynomials, Arblib

setprecision(256)
Arblib.setprecision(256)

lambda = BigFloat("1064e-9")
mu = 4 * big(pi) * BigFloat("1e-7")
epsilon = BigFloat(885418782.0) * BigFloat("1e-20")
c = BigFloat(1 / sqrt(mu * epsilon))
k = 2 * big(pi) / lambda
omega = k * c
K = 2 * k
Zmax = 20 * BigFloat(10^-6)
z0_wave = 0 #5 * 10^-6
z0 = -15 * 10^-6
Q = 0.9 * k
E0 = 1

pmax = 500
pused = 500
#nmax = 47
nmax = 80

lmax = Int(round((Zmax - z0_wave) * k / pi, RoundDown))
lmin = Int(round((-Zmax - z0_wave) * k / pi, RoundUp))
println("lmin: ", lmin, " lmax: ", lmax)

tam = 1000

function load_B_values(filename::String)
    @load filename B1_vals B2_vals B3_vals B4_vals
    return B1_vals, B2_vals, B3_vals, B4_vals
end

# Loads stored B values
B1_data, B2_data, B3_data, B4_data = load_B_values("./GLMT/Bvalues.jld2")

function F(z)
    term1 = exp(Complex{BigFloat}(0, -1) * Q * (z - z0_wave))
    term2 = exp(- 8 * ((z - z0_wave)/Zmax)^2)
    return term1 * term2
end

Fterms = F.( big(pi) .* range(lmin,lmax) ./ k)

function custom_sinint(x)
    x_arb = Arb(x)
    si = sinint(x_arb)
    return BigFloat(si)
end

sinint_table = [custom_sinint((big(pi) * term + k * z0)) for term in (-pmax+lmin):(pmax+lmax)]
println("computed sinint")

function gTM(n, m)
    term1 = - (Complex{BigFloat}(0, 1))^BigFloat(m) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = BigFloat(0)
    for p in range(-pused,pused)
        term3 = BigFloat(0)
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint_table[l - lmin + p + pmax + 1]
        end
        term2 += (B1_data[n, p + pmax + 1] * (m + 1) + B2_data[n, p + pmax + 1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

function gTE(n, m)
    term1 = (Complex{BigFloat}(0, 1))^BigFloat(m + 1) / 2 * (-1.0)^((m - abs(m))/2) * exp(loggamma(BigFloat(n-m+1)) - loggamma(BigFloat(n + abs(m) + 1))) # uses log to prevent overflow with big factorials
    term2 = 0
    for p in range(-pused,pused)
        term3 = 0
        for l in range(lmin, lmax)
            term3 += Fterms[l - lmin + 1] * sinint_table[l - lmin + p + pmax + 1]
        end
        term2 += (B3_data[n, p + pmax + 1] * (m + 1) - B4_data[n, p + pmax + 1] * (1 - m))/2 * term3
    end
    return term1 * term2
end

gTMcoef = [(gTM(n, m)) for n in 1:nmax, m in (-1,1)]
gTEcoef = [(gTE(n, m)) for n in 1:nmax, m in (-1,1)]
println("computed BSCs")
println("gTM(1, 1): ", gTMcoef[1, 2])
println("gTM(2, 1): ", gTMcoef[2, 2])
println("gTM(5, 1): ", gTMcoef[5, 2])
println("gTM(10, 1): ", gTMcoef[10, 2])
println("gTM(25, 1): ", gTMcoef[25, 2])

function cn(n)
    return BigFloat(1.0) / (Complex{BigFloat}(0, k)) * (Complex{BigFloat}(0, -1))^BigFloat(n) * BigFloat(2*n+1) / BigFloat(n * (n+1))
end

cterms = cn.(range(1,nmax))

function custom_besselj(n, x)
    n_acb = Acb(n + 0.5)
    x_acb = Acb(x)
    bess = besselj(n_acb, x_acb)
    result = sqrt(big(pi) / (2 * x)) * Complex{BigFloat}(bess)
    return result
end

function ricatti_bessel(n, x)
    return x * custom_besselj(n, x)
end

function dricatti_besselj(n, x)
    return - n * custom_besselj(n, x) + x * custom_besselj(n-1, x)
end

function dPlm(l,m,theta)
    return central_fdm(3,1)(theta1 -> Plm(l, m, cos(theta1)), theta)
end

function ErTM(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 +=  gTMcoef[n, Int((m + 3)/2)] * exp(Complex{BigFloat}(0,1) * m * phi)
        end
        term1 += cterms[n] * n * (n+1) / r * custom_besselj(n, k * r) * Plm(n, 1, cos(theta)) * term2
    end
    return E0 * term1
end

function EthetaTM(r, theta, phi)
    term1 = BigFloat(0)
    for n in range(1, nmax)
        term2 = BigFloat(0)
        for m in (-1,1)
            term2 +=  gTMcoef[n, Int((m + 3)/2)] * exp(Complex{BigFloat}(0,1) * m * phi)
        end
        term1 += cterms[n] * dricatti_besselj(n, k * r) * dPlm(n,1,theta) * term2
    end
    return E0 * term1 / r
end

function EphiTM(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += BigFloat(m) * gTMcoef[n, Int((m + 3)/2)] * exp(Complex{BigFloat}(0,1) * m * phi)
        end
        term1 += cterms[n] * dricatti_besselj(n, k * r) * Plm(n, 1, cos(theta)) / sin(theta) * term2
    end
    return E0 * term1 * Complex{BigFloat}(0,1) / r
end

@inline function ErTE(r, theta, phi)
    return 0
end

function EthetaTE(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += BigFloat(m) * gTEcoef[n, Int((m + 3)/2)] * exp(Complex{BigFloat}(0,1) * m * phi)
        end
        term1 += cterms[n] * custom_besselj(n, k * r) * Plm(n, 1, cos(theta)) / sin(theta) * term2
    end
    return E0 * term1 * k
end

function EphiTE(r, theta, phi)
    term1 = 0
    for n in range(1, nmax)
        term2 = 0
        for m in (-1,1)
            term2 += gTEcoef[n, Int((m + 3)/2)] * exp(Complex{BigFloat}(0,1) * m * phi)
        end
        term1 += cterms[n] * custom_besselj(n, k * r) * dPlm(n,1,theta) * term2
    end
    return E0 * term1 * Complex{BigFloat}(0,1) * k
end

@inline function Er(r, theta, phi)
    return ErTE(r, theta, phi) + ErTM(r, theta, phi)
end

@inline function Etheta(r, theta, phi)
    return EthetaTE(r, theta, phi) + EthetaTM(r, theta, phi)
end

@inline function Ephi(r, theta, phi)
    return EphiTE(r, theta, phi) + EphiTM(r, theta, phi)
end

function Ex(x, y, z)
    r = BigFloat(sqrt(x^2 + y^2 + (z)^2))
    theta = acos(BigFloat((z)/r))
    phi = atan(BigFloat(y/x))
    comp1 = Er(r, theta, phi) * sin(theta) * cos(phi)
    comp2 = Etheta(r, theta, phi)# * cos(theta) * cos(phi)
    comp3 = - Ephi(r, theta, phi) * sin(phi)
    return comp1 + comp2 + comp3
end

function Ey(x, y, z)
    r = BigFloat(sqrt(x^2 + y^2 + (z)^2))
    theta = acos(BigFloat((z)/r))
    phi = atan(BigFloat(y/x))
    comp1 = Er(r, theta, phi) * sin(theta) * sin(phi)
    comp2 = Etheta(r, theta, phi) * cos(theta) * sin(phi)
    comp3 = Ephi(r, theta, phi) * cos(phi)
    return comp1 + comp2 + comp3
end

function Ez(x, y, z)
    r = BigFloat(sqrt(x^2 + y^2 + (z)^2))
    theta = acos(BigFloat(z/r))
    phi = atan(BigFloat(y/x))
    comp1 = Er(r, theta, phi) * cos(theta)
    comp2 = - Etheta(r, theta, phi) * sin(theta)
    return comp1 + comp2
end

function FW(r, z)
    sum_mackinnon = 0 + 0im
    for n in range(lmin, lmax)
        sum_mackinnon += F(- n * pi / k) * sinc(1/pi * sqrt((k * r)^2 + (z * k + pi * n)^2))
    end
    return sum_mackinnon
end

function S(x)
    return 1/K * reduce(+, Fterms .* exp.(1im * 2 * pi * range(lmin,lmax) * x / K))
end

zvalues = LinRange(-Zmax, Zmax, tam)
x=BigFloat(10^-25)

intensityX = Array{BigFloat}(undef, tam)
Threads.@threads for count in 1:tam
    global intensityX[count] = abs2(Ex(x, 0, zvalues[count] + z0))
end
println("Ex calculated")

intensityZ = Array{BigFloat}(undef, tam)
Threads.@threads for count in 1:tam
    global intensityZ[count] = abs2(Ez(x, 0, zvalues[count] + z0))
end
println("Ez calculated")

#z0_wave = 0
#Fterms = F.( big(pi) .* range(lmin,lmax) ./ k)

intensityFW = [abs2(FW(0, z)) for z in zvalues]
spectra = [abs2(S(x)) for x in range(-k,k,tam)]
println("FW calculated")

intensityX = intensityX / maximum(intensityX)

# plots
fig = Figure(figure_padding=40, size= 3 .*(600, 200), fontsize = 25)
ax1 = Axis(fig[1, 1], title = "(a)",
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.3)
ax2 = Axis(fig[1, 2], title = "(b)",
    xlabel = L"k_z/k", ylabel = L"S(k_z)", aspect=1.3)
ax3 = Axis(fig[1, 3], title = "(c)",
    xlabel = L"z \, (\mu m)", ylabel = L"Magnitude", aspect=1.3)
#ax3 = Axis(fig[1, 3], title = "(c)", xlabel = L"z \, (\mu m)", ylabel = L"ln(erro)", aspect=1.3)

#ax_inset = Axis(fig[1, 1],
#    width=Relative(0.2),
#    height=Relative(0.2),
#    halign=0.05, valign=0.85,
#    xticklabelsvisible=false, yticklabelsvisible=false, yticksvisible=false,
#    ylabelvisible=false, xticksvisible=false)

Fplot = abs2.(F.(zvalues))
zvalues *= 10^6

lines!(ax1, zvalues, Fplot, color="black", label=L"|F(z)|^2")
lines!(ax1, zvalues, intensityFW, color="blue", label=L"|\Psi(0,z)|^2", linestyle=:dash)
lines!(ax2, range(-k,k,tam)/k, spectra, color="black")
#lines!(ax2, zvalues, intensityFW, color="black", label=L"|\Psi(0,z)|^2")
lines!(ax3, zvalues, intensityX, label=L"|E_x|^2", color="blue") #, linestyle=:dash)
lines!(ax3, zvalues, intensityZ, label=L"|E_z|^2", color="red", linestyle=:dashdot)
#lines!(ax3, zvalues, log.(abs.((intensityX .- intensityFW)./intensityX)), color="black")

axislegend(ax1, position=:rt)
#axislegend(ax2, position=:rt)

#colgap!(fig.layout, Relative(0.1))
#rowgap!(fig.layout, 2, 10)
resize_to_layout!(fig)

save("FW_compara.png", fig)
