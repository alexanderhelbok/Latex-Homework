using CairoMakie, GLMakie, Unitful, Roots, LaTeXStrings
import PhysicalConstants.CODATA2018: m_e, ħ

include("init.jl")

CairoMakie.activate!()
GLMakie.activate!()

# load PhysicalConstants
a = 3.81 * u"Å"
V0 = 1 * u"eV"
me = m_e

# equation
J = me*a*V0/ħ^2 |> u"m^(-1)"
K(E) = sqrt(2*me*E)/ħ |> u"m^(-1)"
f(E, k) = J/K(E) * sin(K(E)*a) + cos(K(E)*a) - cos(k*pi)

kvec = range(-1, 1, length=100)

E1 = zeros(length(kvec))
E2 = zeros(length(kvec))

function getE(kvec)
for (i, k) in enumerate(kvec)
    g(E) = f(E, k)
    E1[i] = find_zero(g, 0.1*u"eV")/u"eV"
    E2[i] = find_zero(E -> f(E, k), 7*u"eV")/u"eV"
end
end

E0 = @. uconvert(NoUnits, (ħ^2*(kvec*pi/a)^2/(2me)) / u"eV")

begin
fig = Figure(size=(1200,400), backgroundcolor=:transparent)
ax1 = Axis(fig[1, 1]; 
    kwargs...,
    xlabel = L"k",
    ylabel = L"E\,\, \mathrm{(eV)}",
    xticks = (-2:1:2, [L"-2\pi/a", L"-\pi/a", L"0", L"π/a", L"2\pi/a"]),
    title = "extended scheme",
    titlefont = "Times New Roman",
    titlesize = 24,
    aspect = AxisAspect(1.5),
    backgroundcolor = :transparent
    )
    lines!(ax1, kvec.*2, E0.*4 .+ 0.75, color=:black)

    getE(kvec)
    lines!(ax1, kvec, E1)
    getE(kvec./2 .- 1.5)
    lines!(ax1, kvec./2 .- 1.5, E2, color=:orange)
    getE(kvec./2 .+ 1.5)
    lines!(ax1, kvec./2 .+ 1.5, E2, color=:orange)
    # vertical lines
    vlines!(ax1, [-1,1], color=:gray, linestyle=:dash)

ax2 = Axis(fig[1, 2]; 
    kwargs...,
    xlabel = L"k",
    # ylabel = L"E\ \mathrm{(eV)}",
    xticks = (-1:1:1, [L"-\pi/a", L"0", L"π/a"]),
    title = "reduced scheme",
    titlefont = "Times New Roman",
    titlesize = 24,
    aspect = AxisAspect(1.5),
    backgroundcolor = :transparent
    )
lines!(ax2, kvec, E0 .+ 0.75, color=:black)

getE(kvec)
lines!(ax2, kvec, E1)
lines!(ax2, kvec, E2)

linkyaxes!(ax1, ax2)

# save("V0.pdf", fig)
end
fig