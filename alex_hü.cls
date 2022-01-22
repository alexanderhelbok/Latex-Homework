\author{Alexander Helbok}
\date{\today}
\ProvidesClass{alex_hü}[2021/04/10 Alexander HÜ Class]
\LoadClass[12pt,a4paper]{scrartcl}
\RequirePackage{fullpage}
\RequirePackage[top=2cm, bottom=4.5cm, left=2.5cm, right=2.5cm]{geometry}
%\RequirePackage[utf8]{inputenc} % uncomment if xelatex is not used
\RequirePackage[T1]{fontenc}
\RequirePackage{fontspec}
\RequirePackage[ngerman]{babel}
\RequirePackage{mtpro2}
\RequirePackage{mathtools,amsthm,amsfonts,amscd,siunitx,physics}
\RequirePackage{lastpage}
\RequirePackage{enumerate}
\RequirePackage[headsepline]{scrlayer-scrpage}
\RequirePackage{mathrsfs}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{listings}
\RequirePackage{hyperref}
\RequirePackage{enumitem}
\RequirePackage{multicol}
\RequirePackage{custom}
\RequirePackage{pgfplots}
\RequirePackage{tikz}
\usetikzlibrary{arrows.meta,positioning, calc, shapes, math, decorations.pathreplacing, decorations.markings}
\usepgfplotslibrary{fillbetween}
\pgfplotsset{
	compat=newest,
	legend image code/.code={
		\draw[mark repeat=2,mark phase=2]
		plot coordinates {
			(0cm,0cm)
			(0.25cm,0cm)        %% default is (0.3cm,0cm)
			(0.4cm,0cm)         %% default is (0.6cm,0cm)
		};%
	}
}

%\def\tikzmark#1{\tikz[overlay, remember picture] \node (#1) {};}

%\tikzexternalize[prefix=figures/]

\graphicspath{ {./Graphics/} }

\setmainfont{Times New Roman}

\sisetup{per-mode = symbol, sticky-per}%
\DeclareSIUnit[]\v{\m\per\s}
\DeclareSIUnit[]\a{\m\per\square\s}

\def\vector#1{\hspace{0.1cm} \begin{pmatrix}#1\end{pmatrix}}
\def\dl#1{\underline{\underline{#1}}}	%dl = double underline
%\renewcommand{\vec}{\overrightarrow}

%\renewcommand{\tfrac}[2]{\tfrac{#1}{#2}\hspace{1pt}}
\let\oldsi\si
\renewcommand{\si}[1]{\hspace{4pt}\oldsi{#1}}
\let\oldunit\unit
\renewcommand{\unit}[1]{\hspace{4pt}\oldunit{#1}}

\hypersetup{%
	colorlinks=true,
	linkcolor=blue,
	linkbordercolor={0 0 1}
}

\setlength{\parindent}{0.0in}
\setlength{\parskip}{0.05in}

% define customizeable variables
\newcommand*{\name}[1]{\def\@name{#1}}
\newcommand*{\hwnumber}[1]{\def\@hwnumber{#1}}
\newcommand*{\course}[1]{\def\@course{#1}}


% create header
\clearpairofpagestyles
\setkomafont{pageheadfoot}{\normalfont}
\ihead{\@name}
\chead{\textbf{\Large Blatt \@hwnumber}}
\ohead{\@course \\ \today}
\ofoot{\pagemark}
\headheight 35pt
\headsep 2.5em
\footheight 45pt
