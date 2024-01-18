using CairoMakie, GLMakie, LaTeXStrings

include("init.jl")

CairoMakie.activate!()
GLMakie.activate!()

# load PhysicalConstants
ϵ1 = 0
γ1 = 1
a = 1

kx = range(-2*pi/a, 2*pi/a, length=100)
ky = range(-2*pi/a, 2*pi/a, length=100)

E(kx, ky) = @. ϵ1 + γ1*sqrt(1 + 4*cos(sqrt(3)*ky*a/2)*cos(3*kx*a/2) + 4*cos(3*kx*a/2)^2)

z = [E(x, y) for x in kx, y in ky]

# plot using makie surface
begin
fig = Figure(size=(1200,400))
ax1 = Axis3(fig[1, 1]; 
    # kwargs...,
    xlabel = L"kx",
    ylabel = L"ky",
    zlabel = L"E(k)",
    xticks = ([-2pi,-pi,0,pi,2pi], [L"-2\pi/a", L"-\pi/a", L"0", L"π/a", L"2\pi/a"]),
    yticks = ([-2pi,-pi,0,pi,2pi], [L"-2\pi/a", L"-\pi/a", L"0", L"π/a", L"2\pi/a"]),
    # title = "extended scheme",
    # titlefont = "Times New Roman",
    # titlesize = 24,
    # aspect = AxisAspect(1.5),
    backgroundcolor = :transparent
    )
surface!(ax1, kx, ky, z)

save("Ek.pdf", fig)
fig
end