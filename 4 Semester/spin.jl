using LsqFit, GLMakie, Makie



best, fit = de()

function chisq(f, x, y, p)
    return sum((y - f(x, p)).^2)
end

# f(x, p) = @. p[1]*x + p[2]
# f(x, p) = @. p[1]*x.^2 + p[2]*x + p[3]
f(x, p) = @. p[1]*sin(p[2]*x)

x = [1:0.5:10;]
y = f(x, [1.1, 2.3])
# y = sin.(x) + 0.1*randn(10)

popt = curve_fit(f, x, y, [1.0, 2.0]).param

begin
res1, res2 = 100, 100

a = LinRange(-5, 5, res1)
b = LinRange(-5, 5, res2)

Z = zeros(res1, res2)

@time for i in 1:res1
    for j in 1:res2
        Z[i, j] = chisq(f, x, y, [a[i], b[j]])
    end
end

# find minimum
Min = findmin(Z)[2]
end

println("Minimum at: ", a[Min[1]], ", ", b[Min[2]])
println("Popt at: ", popt[1], ", ", popt[2])

# plot data and fit
begin
fig, ax, s1 = scatter(x, y, color=:blue)
l1 = lines!(ax, 1:0.1:10, f(1:0.1:10, [a[Min[1]], b[Min[2]]]), color=:green)
l2 = lines!(ax, 1:0.1:10, f(1:0.1:10, popt), color=:red)
Legend(fig[1, 2], [l1, l2], ["minimum", "popt"])
fig
end

# plot as wireframe
begin
# fig, ax, w1 = wireframe(a, b, Z, axis=(type=Axis3,), color=:black)
fig, ax, w1 = surface(a, b, Z, axis=(type=Axis3,), color=:black)
scatter!([a[Min[1]]], [b[Min[2]]], [Z[Min[1], Min[2]]], color=:red)
scatter!([popt[1]], [popt[2]], [chisq(f, x, y, popt)], color=:green)
fig
end
