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
\RequirePackage{xifthen}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{ifxetex}
\ifxetex
% Code to be executed only if XeLaTeX is used
\RequirePackage[mtpscr,mtpccal]{mtpro2}
\pdfmapfile{=mtpro2.map}
\fi
\RequirePackage{mathtools, amsthm, amsfonts, siunitx, physics, chemformula, empheq}
\RequirePackage{enumitem, multicol, tcolorbox, array}
\tcbuselibrary{breakable,skins}
\RequirePackage[headsepline]{scrlayer-scrpage}
\RequirePackage{graphicx}
\RequirePackage{listings} %env ignores LATEX commands, keeping all the line breaks (for printing code)
\RequirePackage{hyperref}
\RequirePackage{bookmark, booktabs}
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
\renewcommand*\dd{\mathop{}\!\mathrm{d}}
\renewcommand*{\vector}[1]{\hspace{0.1cm} \begin{pmatrix}#1\end{pmatrix}}
\newcommand*{\iu}{{i\mkern1mu}}
\newcommand*{\dl}[1]{\underline{\underline{#1}}}	%dl = double underline
\newcommand\supfrac[3][]{\mkern-2mu#1\frac{#2}{#3}\rule[-6pt]{0pt}{0pt}}

\newcolumntype{L}{>{$}c<{$}} % math-mode version of "l" column type

% expo function \expo[sign][nominator/denumerator][normal] (all optional)
\ExplSyntaxOn
\NewDocumentCommand{\expo}{>{\SplitArgument{1}{/}} O{} >{\SplitArgument{1}{/}} O{} O{} }{%
	\ifthenelse{\equal{\use_i:nn#1}{-}}
	{\regex_match:nnTF {NoValue}{\use_ii:n#2}
	{\,\mathrm{e}^{\mkern-2mu\use_i:nn#1\rule[-6pt]{0pt}{0pt}\use_i:nn#2}}
	{\,\mathrm{e}^{\mkern-2mu\use_i:nn#1\frac{\use_i:nn#2}{\use_ii:nn#2}\rule[-6pt]{0pt}{0pt}#3}}
	}
	{\regex_match:nnTF {NoValue}{\use_ii:n#1}
	{\,\mathrm{e}^{\mkern-2mu\rule[-6pt]{0pt}{0pt}\use_i:nn#1}}
	{\,\mathrm{e}^{\mkern-2mu\frac{\use_i:nn#1}{\use_ii:nn#1}\rule[-6pt]{0pt}{0pt}\use_i:nn#2}}
	}}
\ExplSyntaxOff

% uint function \uint[Interval]{function}{differential} (interval optional)
\ExplSyntaxOn
\NewDocumentCommand{\uint}{>{\SplitArgument{1}{,}}D[]{-,-} m m}{%
	\ifthenelse{\equal{\use_i:nn#1}{-}}
	{\int#2\, \dd#3}							% case 1: indefinite integral
	{\regex_match:nnTF {NoValue}{\use_ii:n#1}
	{\ifthenelse{\equal{\use_i:nn#1}{o}}
	{\oint #2\, \dd#3}					% case 2: Closed loop Integral
	{\int\c_math_subscript_token{\use_i:nn#1}#2 \dd#3}		% case 3: Volume/Surface Integral
	}
	{\int\limits\c_math_subscript_token{\use_i:nn#1}\c_math_superscript_token{\use_ii:nn#1}#2 \dd#3}	% case 4: Integral with boundaries
	}
}
\ExplSyntaxOff

% redefine cmd to add spacing in math mode
\let\oldtfrac\tfrac
\renewcommand{\tfrac}[2]{\hspace{1pt}\oldtfrac{#1}{#2}\hspace{1pt}}
\let\oldsi\si
\renewcommand{\si}[1]{\hspace{4pt}\oldsi{#1}}
\let\oldunit\unit
\renewcommand{\unit}[1]{\hspace{4pt}\oldunit{#1}}

% --- Override \vec with an invocation of \xvec.
\let\oldvec\vec
\renewcommand{\vec}[1]{\mathbold{#1}}
%\renewcommand{\vec}[1]{\xvec[]{#1}}
% --- Define \dvec and \ddvec for dotted and double-dotted vectors.
\newcommand{\dvec}[1]{\xvec[.]{#1}}
\newcommand{\ddvec}[1]{\xvec[:]{#1}}

\hypersetup{%
	colorlinks=true,
	linkcolor=blue,
	linkbordercolor={0 0 1}
}

% define tcolorbox colors (Latex format very inconsistent)
\definecolor{Title}{HTML}{BFBF00}
\definecolor{Frame}{HTML}{808000}
\definecolor{Background}{HTML}{FFFFF2}

% define tcolorbox env 
\newtcolorbox{mybox}[1]{breakable,pad at break*=1mm, colback=Background,colframe=Frame, fonttitle=\LARGE, colbacktitle=Title, bottomtitle=0.2cm, top=0.35cm, bottom=0.35cm, toptitle=0.2cm, after skip=0.25cm, segmentation style={thick, darkgray!60}, adjusted title={#1}}

\setlength{\parindent}{0.0in}
\setlength{\parskip}{0.05in}

\setlength{\jot}{2ex}

% define customizeable variables
\newcommand*{\name}[1]{\newcommand*\Name{#1}}
\newcommand*{\hwnumber}[1]{\newcommand*\Hwnumber{#1}}
\newcommand*{\course}[1]{\newcommand*\Course{#1}}


% create header
\clearpairofpagestyles
\setkomafont{pageheadfoot}{\normalfont}
\ihead{\Name}
\chead{\textbf{\Large Blatt \Hwnumber}}
\ohead{\Course \\ \today}
\cfoot{\pagemark}
\headheight 35pt
\headsep 2.5em
\footheight 45pt
