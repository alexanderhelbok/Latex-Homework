include("/home/taco/Documents/gp2/Source.jl")
using LaTeXStrings, PhysicalConstants.CODATA2018

import PhysicalConstants.CODATA2018:  ε_0, m_e, e, c_0, N_A

@py import matplotlib.pyplot as plt
@py import matplotlib as mpl
plt.style.use("/home/taco/Documents/gp2/Source.mplstyle")

# mpl.use("pgf")
# plt.rcParams["pgf.texsystem"] = "xelatex"
plt.rc("text", usetex=true)  # enable use of LaTeX in matplotlib
plt.rc("font", family="sans-serif", serif="Times New Roman", size=14)  # font settings
plt.rc("text.latex", preamble="\\usepackage{mtpro2} \\usepackage{siunitx} \\usepackage{physics}")

ρ = 2.7 * u"g/cm^3"
Zal = 13
I(Z) = 16*Z^0.9 * u"eV"
Aal = 26.98 * u"g/mol"
mα = 3727.4 * u"MeV/c^2"
mp = 938.3 * u"MeV/c^2"
mD = 1875.6 * u"MeV/c^2"
Zα = 2
Zp = 1
ZD = 1

T = 2 * u"GeV"
d = 1 * u"cm"

gamma(m, T) = T/(m*c_0^2) + 1
beta(m, T) = sqrt(1 - 1/gamma(m, T)^2)

k = 1/(4π*ε_0)^2 * (4π*e^4/(m_e*c_0^2)) * N_A * ρ * Zal/Aal

dE(m, T, Z) = k * (log(2 * m_e * c_0^2 * beta(m, T)^2 * gamma(m, T)^2/I(Z)) - beta(m, T)^2)

dEα = dE(mα, T, Zα) * d |> u"MeV"
dEp = dE(mp, T, Zp) * d |> u"MeV"
dED = dE(mD, T, ZD) * d |> u"MeV"

Tx = [1:1000:1e6;] * u"MeV"

bgcolor = "#FFFFF2"
# loglogplot
begin
(ax1, ax2) = subplots(2, 1, figsize=(6, 5), sharex=true)[1]
ax1.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mα, Tx, Zα))), label=L"\alpha", color="C0")
ax1.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mp, Tx, Zp))), label=L"p", color="C1")
ax1.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mD, Tx, ZD))), label=L"d", color="C2")

ax2.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mα, Tx, Zα)) / dEα), label=L"\alpha", color="C0")
ax2.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mp, Tx, Zp)) / dEp), label=L"p", color="C1")
ax2.plot(ustrip.(Tx), ustrip.(uconvert.(u"MeV/cm", dE.(mD, Tx, ZD)) / dED), label=L"d", color="C2")
xlabel(L"$T\ (\mathrm{MeV})$")
ax1.set_ylabel(L"-\frac{E}{x}\ (\mathrm{MeV}/\mathrm{cm})")
ax2.set_ylabel(L"-\frac{1}{E} \frac{E}{x}\ (\mathrm{1/cm})")
# yscale("log")
xscale("log")
legend(facecolor=bgcolor)
tight_layout()
# savefig("Graphics/8_loglog.pgf", transparent=true, bbox_inches="tight")
end

# ax = subplots(2, 1, figsize=(6, 3), sharex=true)[0]
# ax[0]
