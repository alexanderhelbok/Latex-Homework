\author{Alexander Helbok}
\date{\today}
\ProvidesClass{alex_summary}[2024/03/24 Alexander Summary Class]
\LoadClass[11pt,a4paper]{scrartcl}
\RequirePackage{fullpage, lastpage}
\RequirePackage[onehalfspacing]{setspace} 
\RequirePackage[left=2.5cm, right=2.5cm, bottom=3.5cm, top=2cm]{geometry}
%\RequirePackage[utf8]{inputenc} % uncomment if xelatex is not used
\RequirePackage[T1]{fontenc}
\RequirePackage{fontspec}
%\setmainfont[Mapping=tex-text]{Times New Roman}
\RequirePackage[backend=biber,
				abbreviate=true,
				doi=true,
				isbn=false,
				url=false,
				style=numeric-comp,
				giveninits=true,
				hyperref=true,
				sorting=none,
				date=year]{biblatex}
\RequirePackage[english]{babel}
\RequirePackage{xifthen, xstring, csquotes, etoolbox}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage[export]{adjustbox}
%\RequirePackage{mtpro2}
%\pdfmapfile{=mtpro2.map}
\RequirePackage{mathtools, amsthm, amsfonts, siunitx, physics, empheq, floatrow, bm}
\RequirePackage{enumitem, multicol, tcolorbox, array, booktabs, vwcol}
\RequirePackage[headsepline]{scrlayer-scrpage}
\RequirePackage[bf, format=plain, font=normalsize, justification=raggedright, singlelinecheck=false]{caption}
\RequirePackage[labelformat=simple]{subcaption}
\RequirePackage[
	pdftitle={},
	pdfsubject={},
	pdfauthor={Alexander Helbok},
	pdfkeywords={}
]{hyperref}
\RequirePackage{graphicx, float, url, wrapfig}
\RequirePackage[nameinlink]{cleveref}
\RequirePackage{bookmark, makecell}
%\RequirePackage{custom}
\RequirePackage{tikz, pgfplots, circuitikz}

\captionsetup[sub]{font=large,labelfont={sf}}
\renewcommand\thesubfigure{(\alph{subfigure})}
\usetikzlibrary{arrows.meta,positioning, calc, shapes, math, decorations.pathreplacing, decorations.markings}
\usepgfplotslibrary{fillbetween}
\tcbuselibrary{breakable,skins}
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

% only print arxiv if no doi is given
\DeclareSourcemap{
	\maps[datatype=bibtex]{
		\map[overwrite]{
			\step[fieldsource=doi, final]
			\step[fieldset=url, null]
			\step[fieldset=eprint, null]
		}  
	}
}

\AtEveryBibitem{\clearfield{pages}}

\Crefname{equation}{Equation}{Equations}
\Crefname{figure}{Figure}{Figures}

\RedeclareSectionCommand[
beforeskip=-\baselineskip,
afterskip=.5\baselineskip]{section}
\makeatletter
\patchcmd{\scr@startchapter}{\if@openright\cleardoublepage\else\clearpage\fi}{}{}{}
\makeatother

%\tikzexternalize[prefix=figures/]

\graphicspath{ {./bilder/} }
\addbibresource{Literatur.bib}

% add custom cite command for thesis drafting
\NewDocumentCommand{\mycite}{m}{
	\colorbox{blue!30}{\citetitle{#1} \autocite{#1}}}

\sisetup{per-mode = symbol, sticky-per}%

% define häufig used cmd
%\renewcommand*\dd{\mathop{}\!\mathrm{d}}
\newcommand*{\iu}{{i\mkern-1mu}}
\newcommand\supfrac[3][]{\mkern-2mu#1\frac{#2}{#3}\rule[-6pt]{0pt}{0pt}}

\newcolumntype{L}{>{$}c<{$}} % math-mode version of "l" column type

% expo function \expo[sign][nominator/denumerator][normal] (all optional)
%BEGIN_FOLD
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

% uint function \uint[Interval]{function}{differential} (interval optional)
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
%END_FOLD

% redefine cmd to add spacing in math mode
\let\oldtfrac\tfrac
\renewcommand{\tfrac}[2]{\hspace{1pt}\oldtfrac{#1}{#2}\hspace{1pt}}
\let\oldunit\unit
\renewcommand{\unit}[1]{\hspace{4pt}\oldunit{#1}}

% define colors
% BEGIN_FOLD
\definecolor{MyBlue}{HTML}{0072B2}
\definecolor{MyOrange}{HTML}{D55E00}
\definecolor{MyRed}{HTML}{F00F0F}
\definecolor{MyGreen}{HTML}{2ca02c}
\definecolor{MyPurple}{HTML}{9467bd}
% END_FOLD

\hypersetup{%
	colorlinks=true,
	breaklinks=true,
	linkcolor=MyBlue,
	urlcolor=MyBlue,
	citecolor=MyBlue,
	linkbordercolor={0 0 1}
}

% place eqn number in brackets for autoref
\makeatletter
\def\tagform@#1{\maketag@@@{\ignorespaces#1\unskip\@@italiccorr}}
\let\orgtheequation\theequation
\def\theequation{(\orgtheequation)}
\makeatother

% create header
% BEGIN_FOLD
\clearpairofpagestyles
\setkomafont{pageheadfoot}{\normalfont}
\ihead{Alexander Helbok}
\chead{}
\ohead{\today}
\headheight \baselineskip
\headsep 2.5em
\footheight 45pt
% END_FOLD

\newtcolorbox{question1}{
	enhanced,
	adjusted title=\vspace{-1cm}Questions,
	fonttitle=\bfseries,
	halign title=center,
	colbacktitle=white,
	coltitle=black,
	titlerule=0pt,
	sharp corners=all,
	colback=white,
	colframe=orange,
	toprule=0pt,
	bottomrule=0pt,
	leftrule=1.5pt,
	rightrule=0pt,
	overlay={\draw[orange,line width=1.5pt] (frame.north west) -- ++(1cm,0pt);}
}

\newtcolorbox{question2}{
	enhanced,
	sharp corners=all,
	colback=white,
	colframe=orange,
	toprule=0pt,
	bottomrule=0pt,
	leftrule=1.5pt,
	rightrule=1.5pt
}

\newtcolorbox{question3}{
	enhanced,
	sharp corners=all,
	colback=white,
	colframe=orange,
	toprule=0pt,
	bottomrule=0pt,
	leftrule=0pt,
	rightrule=1.5pt,
	overlay={\draw[orange,line width=1.5pt] (frame.south east) -- ++(-1cm,0pt);}
}