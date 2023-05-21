include("/home/taco/Documents/gp2/Source.jl")
using LaTeXStrings

@py import matplotlib.pyplot as plt
@py import matplotlib as mpl
plt.style.use("/home/taco/Documents/gp2/Source.mplstyle")

# mpl.use("pgf")
# plt.rcParams["pgf.texsystem"] = "xelatex"
plt.rc("text", usetex=true)  # enable use of LaTeX in matplotlib
plt.rc("font", family="sans-serif", serif="Times New Roman", size=14)  # font settings
plt.rc("text.latex", preamble="\\usepackage{mtpro2} \\usepackage{siunitx}")

p = [0.155, 0.250, 0.350, 0.415, 0.450, 0.483, 0.517, 0.533, 0.550, 0.567, 0.583] * u"bar"
R = [36.2, 35.5, 36.0, 35.3, 35.9, 34.5, 30.4, 16.6, 10.5, 4.8, 1.7] * u"s^-1"

x0 = 6 * u"cm"

x = p * x0 * u"bar^-1"

model(x, p) = @. p[2] * (1 + exp(-2*p[1]*(x - p[4])))^(-1) + p[3]

popt, pcov, ci = bootstrap(model, ustrip.(x[1:end]), ustrip.(R[1:end]), p0=[2.5, -60.0, 35.0, 3.4])

popt = curve_fit(model, ustrip.(x[1:end]), ustrip.(R[1:end]), [2.5, -60.0, 35.0, 3.4]).param

point = -log((-10*popt[2])/(9*popt[3]) - 1)/(2*popt[1]) + popt[4]

begin
figure(figsize=(6, 4))
scatter(x, R, label=L"$\mathrm{Messwerte}$")
scatter(point, model(point, popt), facecolor="none", edgecolor="k", zorder=5, s=80, lw=2)
plot(ci.x, model(ci.x, popt), color="red", label=L"$f(x) \propto (1 + \mathrm{e}^{-2kx})^{-1} $")
axvline(point, color="k", ls="--", zorder=0)
axhline(model(point, popt), color="k", ls="--",  zorder=0)
xlabel(L"$\mathrm{x_{eff}\ (cm)}$")
ylabel(L"$\mathrm{R\ (1/s)}$")
# grid()
legend()
tight_layout()
# savefig("plot1.pdf", transparent=true, bbox_inches="tight")
end


Ekin = (point/0.31)^(2/3) * u"MeV"

model2(x, p) = @. (2 * p[1] * p[2] * exp(-2*p[1]*(x - p[4])))/(1 + exp(-2*p[1]*(x - p[4])))^2

# also use numerical differentiation
x2 = ustrip.(x[1:end-1] + diff(x)/2)
dR = ustrip.(diff(R))./ustrip.(diff(x))

# Zeichnen Sie schließlich ein Diagramm der Absorptionsrate (also der Anderung der gemessenen Rate) der α-Teilchen als Funktion der L¨ange x bei Normaldruck. Interpretieren Sie das Ergebnis.
begin
scatter(x2, dR)
plot(ci.x, model2(ci.x, popt))
xlabel("x (cm)")
ylabel(L"\frac{dR}{dx} (s^{-1} cm^{-1})")
tight_layout()
# savefig("dRdx.pdf", transparent=true, bbox_inches="tight")
end