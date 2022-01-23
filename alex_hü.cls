\author{Alexander Helbok}
\date{\today}
\ProvidesClass{alex_hü}[2021/04/10 Alexander HÜ Class]
\LoadClass[12pt,a4paper]{scrartcl}
\RequirePackage{fullpage, lastpage}
\RequirePackage[top=2cm, bottom=4.5cm, left=2.25cm, right=2.25cm]{geometry}
%\RequirePackage[utf8]{inputenc} % uncomment if xelatex is not used
\RequirePackage[T1]{fontenc}
\RequirePackage{fontspec}
\RequirePackage[ngerman]{babel}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{mtpro2}
\RequirePackage{mathtools,amsthm,amsfonts,siunitx,physics}
\RequirePackage{enumitem, multicol, tcolorbox}
\RequirePackage[framemethod=tikz]{mdframed} %for \Sepline cmd in tcolorbox
\tcbuselibrary{breakable}
\RequirePackage[headsepline]{scrlayer-scrpage}
\RequirePackage{graphicx}
\RequirePackage{listings} %env ignores LATEX commands, keeping all the line breaks (for printing code)
\RequirePackage{hyperref}
\RequirePackage{custom}
\RequirePackage{tikz, pgfplots}
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

%\tikzexternalize[prefix=figures/]

\graphicspath{ {./Graphics/} }

\setmainfont{Times New Roman}

\sisetup{per-mode = symbol, sticky-per}%
\DeclareSIUnit[]\v{\m\per\s}
\DeclareSIUnit[]\a{\m\per\square\s}

% define häufig used cmd
\renewcommand{\vector}[1]{\hspace{0.1cm} \begin{pmatrix}#1\end{pmatrix}}
\newcommand{\dl}[1]{\underline{\underline{#1}}}	%dl = double underline

% redefine cmd to add spacing in math mode
\let\oldtfrac\tfrac
\renewcommand{\tfrac}[2]{\hspace{1pt}\oldtfrac{#1}{#2}\hspace{1pt}}
\let\oldsi\si
\renewcommand{\si}[1]{\hspace{4pt}\oldsi{#1}}
\let\oldunit\unit
\renewcommand{\unit}[1]{\hspace{4pt}\oldunit{#1}}

\hypersetup{%
	colorlinks=true,
	linkcolor=blue,
	linkbordercolor={0 0 1}
}

% define \Sepline cmd for use in tcolorbox
\gdef\Sepline{%
	\par\noindent\makebox[\linewidth][l]{%
		\hspace*{-\mdflength{innerleftmargin}-1.24cm}%
		\tikz\draw[thick,dashed,darkgray!60] (0,0) --%
		(\textwidth+\the\mdflength{innerleftmargin}+\the\mdflength{innerrightmargin}+0.3cm,0);
	}\par\nobreak}

% define tcolorbox env 
\newtcolorbox{mybox}[1]{breakable,pad at break*=1mm, colback=yellow!5!white,colframe=yellow!50!black, fonttitle=\LARGE, colbacktitle=yellow!75!black, bottomtitle=0.2cm, top=0.35cm, bottom=0.35cm, toptitle=0.2cm, adjusted title={#1}}

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
