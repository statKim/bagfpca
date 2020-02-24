%%==========================================================================
%% Style File for Communications for Statistical Applications and Methods (CSAM)
%% Prepared by the Editorial Board of CKSS (2008-2009)
%% Korean Statistical Society (office@kss.or.kr)
%%==========================================================================
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{ckss}
              [2008/11/12 v0.1 Communications for Statistical Applications and Methods class]            
\RequirePackage{ifthen}
\newif\if@korean            
%%=============
%% Page layout
%%=============
\DeclareOption{a4paper}
   {\setlength\paperheight    {297mm}%
    \setlength\paperwidth     {210mm}%
    \setlength\textheight     {19.5cm}%
    \setlength\textwidth      {14.2cm}%
    \setlength\oddsidemargin  {0.90 cm}%
    \setlength\evensidemargin {0.90 cm}%
    \setlength\parindent      {1.5em}%
    \setlength\headheight     {20 pt}%
    \setlength\textfloatsep{8pt plus 2pt minus 2pt}%
    \setlength\textfloatsep{8pt plus 2pt minus 2pt}}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
\DeclareOption{kor}{\@koreantrue}
\DeclareOption{eng}{\@koreanfalse}
\input{size10.clo}
\ExecuteOptions{a4paper,kor,final}
\ProcessOptions
\RequirePackage{amsmath}
\RequirePackage{txfonts}%{mathptm}
\RequirePackage{fancyhdr}
\RequirePackage{amssymb}       
\RequirePackage[format=plain,indention=0cm,labelsep=colon,%
                singlelinecheck=false,font=small,labelfont={bf,sf},%
                textfont={sl},parskip=4pt]{caption}
\captionsetup[figure]{skip=2pt,labelfont={sf},textfont={it},justification=centerlast}
\captionsetup[table]{skip=2pt,labelfont={sf},textfont={rm},justification=raggedright}
\if@korean
   \RequirePackage[nojosa]{hangul}
   \def\ksnamedef#1#2{%
   \AtBeginDocument{\@namedef{ks#1}{#2}
                    \@namedef{#1}{\expandafter\protect\csname ks#1\endcsname}}}
        \newcommand{\nohyphens}{\hyphenpenalty=10000\exhyphenpenalty=10000\relax}
        \def\appendixname{부록}
  \ksnamedef{tablename}{Table}
  \ksnamedef{abstractname}{Abstract}        
  \ksnamedef{figurename}{Figure}        
  \ksnamedef{refname}{References}
        \def\theoremname{\textgt{정리}}
        \def\corollaryname{\textgt{보조정리}}
        \def\definitionname{\textgt{정의}}
        \def\proofname{\textbf{증명}}
        \def\propositionname{\textgt{명제}}
        \def\assumptionname{\textrm{가정}}
        \def\examplename{\textrm{예제}}
        \def\keywordsname{Keywords:}
        \def\@journaltitle{\makebox[80pt][s]{\fontsize{7}{7}\selectfont%
        \textgt{한 국 통 계 학 회 논 문 집}}}
        \newcommand{\sff}{\그래프}
        \newcommand{\thmshape}{\rmfamily}
        \def\@kesect{\sff\bfseries\fontsize{10}{10}\selectfont}
        \def\@kesubsect{\sff\fontsize{10}{10}\selectfont}
        \def\@kesubsubsect{\itshape\fontsize{9.5}{10}\selectfont}
        \def\@ixpt{9}
        \def\@xpt{9.5}
        \def\@xipt{10.45}
        \renewcommand\normalsize{%
           \@setfontsize\normalsize\@xpt\@xiipt
           \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
           \abovedisplayshortskip \z@ \@plus3\p@
           \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
           \belowdisplayskip \abovedisplayskip
           \let\@listi\@listI}
        \normalsize
        \renewcommand\small{%
        \@setfontsize\small\@ixpt{11}%
        \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus2\p@
        \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
        \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
        \belowdisplayskip \abovedisplayskip}
        \renewcommand\footnotesize{%
              \@setfontsize\footnotesize\@viiipt{9.5}%
              \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
              \abovedisplayshortskip \z@ \@plus\p@
              \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
              \def\@listi{\leftmargin\leftmargini
                   \topsep 3\p@ \@plus\p@ \@minus\p@
                   \parsep 2\p@ \@plus\p@ \@minus\p@
                   \itemsep \parsep}%
              \belowdisplayskip \abovedisplayskip}
        \renewcommand{\baselinestretch}{1.2}              
\else
        \def\appendixname{Appendix}
        \def\tablename{Table}
        \def\figurename{Figure}
        \def\abstractname{Abstract}
        \def\keywordsname{Keywords:}
        \def\refnamex{References}
        \def\@journaltitle{\fontsize{7}{7}\selectfont Communications for %
                          Statistical Applications and Methods}
        \def\theoremname{\fontsize{10}{10}\usefont{T1}{phv}{b}{n}\selectfont Theorem}
        \def\lemmaname{\fontsize{10}{10}\usefont{T1}{phv}{b}{n}\selectfont Lemma}
        \def\corollaryname{\fontsize{10}{10}\usefont{T1}{phv}{b}{n}\selectfont Corollary}
        \def\definitionname{\textbf{Definition}}
        \def\propositionname{\fontsize{10}{10}\usefont{T1}{phv}{b}{n}\selectfont Proposition}
        \def\assumptionname{\sffamily Assumption}
        \def\proofname{\textbf{Proof}}
        \def\examplename{\usefont{T1}{cmr}{bx}{it} \selectfont Example}
        \def\remarkname{\usefont{T1}{cmr}{bx}{it} \selectfont Remark}
        \def\sff{\sffamily}
        \newcommand{\thmshape}{\itshape}
        \renewcommand{\baselinestretch}{1}
        \def\@kesect{\sff\bfseries\fontsize{11}{11}\selectfont}
        \def\@kesubsect{\sff\fontsize{11}{11}\selectfont}
        \def\@kesubsubsect{\slshape\fontsize{10}{10}\selectfont}
        \captionsetup[figure]{skip=0pt, labelfont={sf},textfont={sl},
                             justification=centerlast}
        \captionsetup[table]{skip=0pt, labelfont={sf},textfont={rm},
                            justification=raggedright}
        \def\convert#1{\def\@month{#1}
            \ifnum\@month=1  {January} \fi
            \ifnum\@month=2  {February} \fi
            \ifnum\@month=3  {March} \fi
            \ifnum\@month=4  {April} \fi
            \ifnum\@month=5  {May} \fi
            \ifnum\@month=6  {June} \fi
            \ifnum\@month=7  {July} \fi
            \ifnum\@month=8  {August} \fi
            \ifnum\@month=9  {September} \fi
            \ifnum\@month=10 {October} \fi
            \ifnum\@month=11 {November} \fi
            \ifnum\@month=12 {December} \fi}
        \def\@ixpt{9}
        \def\@xpt{10}
        \def\@xipt{10.95}
        \renewcommand\normalsize{%
           \@setfontsize\normalsize\@xpt\@xiipt
           \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
           \abovedisplayshortskip \z@ \@plus3\p@
           \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
           \belowdisplayskip \abovedisplayskip
           \let\@listi\@listI}
        \normalsize
        \renewcommand\small{%
        \@setfontsize\small\@ixpt{11}%
        \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus2\p@
        \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
        \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
        \belowdisplayskip \abovedisplayskip}
        \renewcommand\footnotesize{%
              \@setfontsize\footnotesize\@viiipt{9.5}%
              \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
              \abovedisplayshortskip \z@ \@plus\p@
              \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
              \def\@listi{\leftmargin\leftmargini
                   \topsep 3\p@ \@plus\p@ \@minus\p@
                   \parsep 2\p@ \@plus\p@ \@minus\p@
                   \itemsep \parsep}%
              \belowdisplayskip \abovedisplayskip}
        \renewcommand{\baselinestretch}{1}              
\fi

\setlength\leftmargini  {1.5em}
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty

\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{3}%2
\renewcommand\topfraction{1}%0.7
\setcounter{bottomnumber}{3}%1
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{5}%3
\renewcommand\textfraction{0}%.2
\renewcommand\floatpagefraction{1}%0.5
%\theoremstyle{break}
\newtheorem{theorem}{\theoremname}%[section]
\newtheorem{corollary}{\corollaryname}%[section]
\newtheorem{definition}{\definitionname}%[section]
\newtheorem{proposition}{\propositionname}%[section]
\newtheorem{assumption}{\assumptionname}%[section]
\newtheorem{lemma}{\lemmaname}%[section]
\setcounter{secnumdepth}{3}
\newcounter{endpage}
\newcounter{remark}
\newcounter{example}
\newcounter{ref}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
%\newcounter {paragraph}[subsubsection]
%\newcounter {subparagraph}[paragraph]
\renewcommand\thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
%\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
%\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\renewcommand{\theref}{}
\renewcommand{\thetheorem}{\@arabic\c@theorem{}}
\renewcommand{\theproposition}{\@arabic\c@proposition{}}
\renewcommand{\thecorollary}{\@arabic\c@corollary{}}
\renewcommand{\thelemma}{\@arabic\c@lemma{}}
\renewcommand{\theassumption}{\@arabic\c@assumption{}}
\renewcommand{\theexample}{\arabic{example}{}}
\renewcommand{\theremark}{\arabic{remark}{}}
\renewcommand{\thedefinition}{\arabic{definition}{}}
%\renewcommand{\thesection}{\arabic{section}}
%\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}%
\DeclareSymbolFont{AMSa}{U}{msa}{m}{n}
\def\m@th{\mathsurround\z@}
\def\hb@xt@{\hbox to}
\DeclareRobustCommand*\textsuperscript[1]{%
  \@textsuperscript{\selectfont#1}}
\def\@textsuperscript#1{%
  $\m@th^{\mbox{\fontsize\sf@size\z@$#1$}}$}

% For sample file
\def\AMSLaTeX{$\cal A\kern-.1667em \lower.5ex\hbox{$\cal M$}%
              \kern-.125em S$-\LaTeX}
% SUPERSCRIPTS AND SUBSCRIPTS
% Declaration which causes superscripts and subscripts to be output at normal
% text size
\long\def\@footnotetext#1{\insert\footins%
 {\def\baselinestretch{1}\footnotesize%
 \interlinepenalty\interfootnotelinepenalty%
 \splittopskip\footnotesep \splitmaxdepth%
 \dp\strutbox \floatingpenalty \@MM \hsize\columnwidth%
 \@parboxrestore \edef%
 \@currentlabel{\csname p@footnote\endcsname\@thefnmark}
 \@makefntext{\rule{\z@}{\footnotesep}\ignorespaces #1\strut}}}

% Singlespaced hanging paragraph environment for uncited bibliographies.
\def\singlespace{\vskip\baselineskip \def\baselinestretch{1}\@normalsize
      \vskip-\baselineskip}
\def\endsinglespace{\par}
\if@korean
     % Indentation after sectioning command
        \let\@afterindentfalse\@afterindenttrue
        \@afterindenttrue
\fi
\def\ps@centerheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{{\footnotesize \thepage}\hfill\leftmark}%
    \def\@oddhead{\rightmark\hfill {\footnotesize\thepage}}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble}
\def\ps@onesideheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\hfil\slshape\fontsize{8}{9}\selectfont\leftmark}%
    \def\@oddhead{\rightmark \hfil}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble}
\def\\numberwithinn#1#2{\@ifundefined{c@#1}{\@nocounterr{#1}}{%
  \@ifundefined{c@#2}{\@nocnterr{#2}}{%
  \@addtoreset{#1}{#2}%
  \toks@\expandafter\expandafter\expandafter{\csname the#1\endcsname}%
  \expandafter\xdef\csname the#1\endcsname
    {\expandafter\noexpand\csname the#2\endcsname
     .\the\toks@}}}}
\def\@begintheorem#1#2{\trivlist
   \item[\hskip \labelsep{\bf #1\ #2\@addpunct{.}}]\ignorespaces \thmshape}%\itshape}
\def\@opargbegintheorem#1#2#3{\trivlist
      \item[\hskip \labelsep{\bfseries #1\ #2\@addpunct{.} (#3)}\ignorespaces]\thmshape}%\itshape}
\def\pageref#1{\expandafter\@setref\csname r@#1\endcsname
                                   \@secondoftwo{#1}}
\newenvironment{proof}[1][\proofname]{
  \trivlist
  \item[\hskip\labelsep \textbf{#1}:]\ignorespaces
}  
%\def\proof[#1]{\trivlist
%           \item[\hskip \labelsep {\reset@font {\textbf{#1} }}]}
           
\def\endproof{{\hspace*{\fill}\nolinebreak[2]\hspace*{\fill}}$\Box$\endtrivlist}
\\numberwithinn{equation}{section}
\def\@seccntformat#1{\csname the#1\endcsname. {}}%\hspace{0.5em}}
\newcommand\section{\@startsection
         {section}{1}{0mm}%
         {-1.2\baselineskip}%
         {0.5\baselineskip}%
         {\raggedright \@kesect}}%\sff\bfseries\fontsize{10}{10}\selectfont}}%
\newcommand\subsection{\@startsection
         {subsection}{2}{0mm}%
         {-\baselineskip}%
         {0.4\baselineskip}%
         {\raggedright \@kesubsect}}%\sff\bfseries\fontsize{9}{9}\selectfont}}%
\newcommand\subsubsection{\@startsection
         {subsubsection}{3}{0mm}%
         {-\baselineskip}%
         {0.4\baselineskip}%
         {\raggedright \@kesubsubsect}}%\sff\bfseries\fontsize{9}{9}\selectfont}}%
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
%\renewcommand \theequation {\@arabic\c@equation}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename\nobreakspace\thefigure}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
%\newlength\abovecaptionskip
%\newlength\belowcaptionskip
%\setlength\abovecaptionskip{10\p@}
%\setlength\belowcaptionskip{0\p@}
%\long\def\@makecaption#1#2{%
%  \vskip\abovecaptionskip
%  \sbox\@tempboxa{#1: #2}%
%  \ifdim \wd\@tempboxa >\hsize
%    #1: #2\par
%  \else
%    \global \@minipagefalse
%    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%  \fi
%  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
     {\section*{\refname}%
      \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\long\def\@makefntext#1{\parindent 1em\noindent
    \hbox to 4pt{\footnotesize{$^{\@thefnmark}$}}\hskip 1pt\parbox[t]{0.95\linewidth}%
                               {\normalfont #1}}%
                               \setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\def\@final{final}
\def\submit#1{\def\@submit{#1}
             \ifx\@submit\@final \@@final
             \else \@@draft \fi}
\def\@@final{\@twosidetrue \@mparswitchtrue
        \def\@DOI{\fontsize{7}{7}\selectfont DOI: http://dx.doi.org/10.5351/CSAM.\@year.\@volumn.\@ho.\dpage}
        \if@korean
          \def\@journalhead{\makebox[80pt][s]{\fontsize{7}{7}\selectfont \@year,%
                            \hss \@volumn권, \hss \@ho호, \hss \thepage--\pageref{sec:\theendpage}\phantom{.}}}
          \def\heading##1##2{\def\@righttext{\fontsize{8}{8}\selectfont\textgt{##1}}
                             \def\@lefttext{\fontsize{8}{8}\selectfont\textgt{##2}}}
        \else
          \def\@journalhead{\fontsize{7}{7}\selectfont \@year, Vol.
                            \@volumn, No. \@ho, \thepage--\pageref{sec:\theendpage}}
          \def\heading##1##2{\def\@righttext{\fontsize{9}{9}\selectfont\textsf{##1}}
                             \def\@lefttext{\fontsize{9}{9}\selectfont\textsf{##2}}}
        \fi        
        \pagestyle{centerheadings}
        \def\receive##1##2##3{\def\@receiveyear{##1}\def\@receivemonth{##2}\def\@receiveday{##3}}
        \def\revise##1##2##3{\def\@reviseyear{##1}\def\@revisemonth{##2}\def\@reviseday{##3}}
        \def\accept##1##2##3{\def\@acceptyear{##1}\def\@acceptmonth{##2}\def\@acceptday{##3}}}
\def\@@draft{\@twosidefalse \@mparswitchfalse
            \pagestyle{centerheadings}
           \def\heading##1{\markboth{}{##1}}}
                      
\def\title#1{\cleardoublepage
    \def\@title{#1}\par
    \renewcommand\thefootnote{{}}%{\@fnsymbol\c@footnote}%
    \setcounter{footnote}{1}
    \refstepcounter{endpage}
    \if@twocolumn
       \ifnum \col@number=\@ne
          \@maketitle
       \else
          \twocolumn[\@maketitle]%
       \fi
    \else
       \newpage
       \global\@topnum\z@   % Prevents figures from going at top of page.
       \@maketitle
    \fi
    \setcounter{footnote}{0}
}
\def\etitle#1{%
    \@koreanfalse
        \def\@ixpt{9}
        \def\@xpt{10}
        \def\@xipt{10.95}
        \renewcommand\normalsize{%
           \@setfontsize\normalsize\@xpt\@xiipt
           \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
           \abovedisplayshortskip \z@ \@plus3\p@
           \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
           \belowdisplayskip \abovedisplayskip
           \let\@listi\@listI}
%        \normalsize
        \renewcommand\small{%
        \@setfontsize\small\@ixpt{11}%
        \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
        \abovedisplayshortskip \z@ \@plus2\p@
        \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
        \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
        \belowdisplayskip \abovedisplayskip}
        \renewcommand\footnotesize{%
              \@setfontsize\footnotesize\@viiipt{9.5}%
              \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
              \abovedisplayshortskip \z@ \@plus\p@
              \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
              \def\@listi{\leftmargin\leftmargini
                   \topsep 3\p@ \@plus\p@ \@minus\p@
                   \parsep 2\p@ \@plus\p@ \@minus\p@
                   \itemsep \parsep}%
              \belowdisplayskip \abovedisplayskip}
        \renewcommand{\baselinestretch}{1}                  
%    \renewcommand{\sff}{\sffamily}
    \def\abstractname{Abstract}
    \def\keywordsname{Keywords:}
    \def\@etitle{#1}\par
    \setcounter{footnote}{1}
    \renewcommand\thefootnote{}%{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\@thefnmark}}}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@makeetitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
     \@makeetitle
    \fi}
\def\@maketitle{%
   \if@twoside
       \thispagestyle{fancy}
       \fancyfoot{}
       \renewcommand{\headrulewidth}{0pt}
       \rhead{\@DOI\\[-0.4em] \phantom{\@journalhead}}
       \if@korean
          \lhead{\@journaltitle \\[-0.4em] \@journalhead}
       \else
          \lhead{\@journaltitle \\[.2em] \@journalhead}
       \fi
       \markboth{\@lefttext}{\@righttext}
   \else
       \thispagestyle{empty}
   \fi
   \null
   \vspace*{-0.2in}%
   \renewcommand{\baselinestretch}{1.4} \@normalsize
   \begin{center}
   \if@korean
       {\fontsize{14}{14}\bf\sff\selectfont \@title} \par%
   \else
       \vspace*{-1.18em}
       {\fontsize{14.5}{14.5}\usefont{T1}{phv}{b}{n}\selectfont \@title} \par
    \fi
    \end{center}
    \vskip 1.2em%%
}
\def\@makeetitle{%
    \def\@evenhead{\footnotesize\thepage\hfil\phantom{\thepage}}%
    \def\@oddhead{\footnotesize\phantom{\thepage}\hfil\thepage}%
  \null
   \vspace*{-0.256in}%
   \renewcommand{\baselinestretch}{1.4} \@normalsize
   \begin{center}
       \vspace*{-1.18em}
       {\fontsize{14.5}{14.5}\usefont{T1}{phv}{b}{n}\selectfont \@etitle} \par
    \end{center}
    \vskip 1.2em%%
}%
\def\author#1{
    \renewcommand{\baselinestretch}{1.2} \@normalsize
    \def\@author{#1}
        \setcounter{footnote}{0}
    \renewcommand{\thefootnote}{\arabic{footnote}}%{\@fnsymbol\c@footnote}
    \def\@makefnmark{\hspace{4pt} \llap{\@textsuperscript{\normalfont\@thefnmark,}}}%
    \begin{center}
       \if@korean
           {\sffamily\fontsize{10.5}{12} \selectfont \@author}
       \else
           {\sffamily\fontsize{11}{12} \selectfont \@author}
       \fi    
    \end{center} 
    \@thanks
    \theendnotes
}
\def\HRule{\noindent\rule{\textwidth}{0.6pt}}
\def\abstract{%
%    \normalsize
    \HRule
    \par 
    \vspace*{-2pt}
%    \if@korean
%        \renewcommand{\baselinestretch}{1.1}
%    \else    
        \renewcommand{\baselinestretch}{1} 
%    \fi    
        \centerline{\sffamily\bfseries \abstractname} %\par \noindent %\vspace{-0.4em}%{\vspace{\z@}}%-1.0
   \vspace{2pt}
	\small
    \setcounter{footnote}{0}
    \renewcommand{\thefootnote}{\arabic{footnote}}
}
\def\endabstract{\par}

%gogo quote to quotation
\def\ekeywords#1{
        \par\vspace{5pt}
        \renewcommand{\baselinestretch}{1}
        \normalsize
        \noindent 
        {\sffamily\keywordsname{} #1}
        \par \vspace*{-5pt}
        \HRule
        \normalfont \normalsize
        \ifx\@submit\@final\label{sec:\theendpage} \fi
        \vfill 
}
\def\keywords#1{
        \par\vspace{5pt}
        \if@korean
           \renewcommand{\baselinestretch}{1.1}
        \else
           \renewcommand{\baselinestretch}{1}
        \fi
           \normalsize
           \noindent
           \sffamily\keywordsname{} #1
        \par \vspace*{-5pt}
        \HRule
        \normalfont\normalsize
        \medskip
}
\newcommand{\Appendix}{\@startsection%
         {section}{1}{0mm}%
         {-\baselineskip}%
         {0.5\baselineskip}%
         {\raggedright \@kesect \appendixname{} \thesection: \hskip 0.3em}}
\newcommand{\sAppendix}{%
            \@startsection {section}{1}{0mm}%\z@}%
             {-\baselineskip}%
             {0.5\baselineskip}%
             {\raggedright \@kesect \appendixname: }}
\newcommand{\appendix}{\par
\\numberwithinn{table}{section}
\\numberwithinn{figure}{section}
  \setcounter{section}{0}%
  \setcounter{equation}{0}
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}
  \setcounter{table}{0}
  \renewcommand{\thefigure}{\@Alph\c@section.\arabic{figure}}
%  \renewcommand{\thefigure}{\Alph{figure}}
  \renewcommand{\thetable}{\@Alph\c@section.\arabic{table}}
  \renewcommand{\thesection}{\@Alph\c@section}
  \renewcommand{\thesubsection}{\@Alph\c@section.\arabic{subsection}}%
  \def\@seccntformat##1{\csname the##1 \endcsname {}}
  \renewcommand{\section}{%
                \secdef\Appendix\sAppendix}}
\def\remark{\stepcounter{remark}
    \trivlist
    \item[\hskip \labelsep {\reset@font \bfseries\itshape \remarkname{} \theremark\@addpunct{.}}]\ignorespaces}
\def\endremark{\endtrivlist}
\def\example{\stepcounter{example}
    \trivlist
    \item[\hskip \labelsep {\reset@font \bfseries\itshape \examplename{} \theexample\@addpunct{.}}]\ignorespaces}
\def\endexample{\endtrivlist}

\def\volumn#1#2#3{\def\@volumn{#1}\def\@ho{#2}\def\@year{#3}}
\def\setpagenum#1{\def\pagenum{#1} \setcounter{page}{\pagenum}
       \def\dpage{\pagenum}
       \ifnum\pagenum=1  \def\dpage{001} \fi \ifnum\pagenum=3  \def\dpage{003} \fi \ifnum\pagenum=5  \def\dpage{005} \fi 
       \ifnum\pagenum=7  \def\dpage{007} \fi \ifnum\pagenum=9  \def\dpage{009} \fi \ifnum\pagenum=11 \def\dpage{011} \fi 
       \ifnum\pagenum=13 \def\dpage{013} \fi \ifnum\pagenum=15 \def\dpage{015} \fi \ifnum\pagenum=17 \def\dpage{017} \fi 
       \ifnum\pagenum=19 \def\dpage{019} \fi \ifnum\pagenum=21 \def\dpage{021} \fi \ifnum\pagenum=23 \def\dpage{023} \fi 
       \ifnum\pagenum=25 \def\dpage{025} \fi \ifnum\pagenum=27 \def\dpage{027} \fi \ifnum\pagenum=29 \def\dpage{029} \fi 
       \ifnum\pagenum=31 \def\dpage{031} \fi \ifnum\pagenum=33 \def\dpage{033} \fi \ifnum\pagenum=35 \def\dpage{035} \fi 
       \ifnum\pagenum=37 \def\dpage{037} \fi \ifnum\pagenum=39 \def\dpage{039} \fi \ifnum\pagenum=41 \def\dpage{041} \fi 
       \ifnum\pagenum=43 \def\dpage{043} \fi \ifnum\pagenum=45 \def\dpage{045} \fi \ifnum\pagenum=47 \def\dpage{047} \fi 
       \ifnum\pagenum=49 \def\dpage{049} \fi \ifnum\pagenum=51 \def\dpage{051} \fi \ifnum\pagenum=53 \def\dpage{053} \fi 
       \ifnum\pagenum=55 \def\dpage{055} \fi \ifnum\pagenum=57 \def\dpage{057} \fi \ifnum\pagenum=59 \def\dpage{059} \fi 
       \ifnum\pagenum=61 \def\dpage{061} \fi \ifnum\pagenum=63 \def\dpage{063} \fi \ifnum\pagenum=65 \def\dpage{065} \fi 
       \ifnum\pagenum=67 \def\dpage{067} \fi \ifnum\pagenum=69 \def\dpage{069} \fi \ifnum\pagenum=71 \def\dpage{071} \fi 
       \ifnum\pagenum=73 \def\dpage{073} \fi \ifnum\pagenum=75 \def\dpage{075} \fi \ifnum\pagenum=77 \def\dpage{077} \fi 
       \ifnum\pagenum=79 \def\dpage{079} \fi \ifnum\pagenum=81 \def\dpage{081} \fi \ifnum\pagenum=83 \def\dpage{083} \fi
       \ifnum\pagenum=85 \def\dpage{085} \fi \ifnum\pagenum=87 \def\dpage{087} \fi \ifnum\pagenum=89 \def\dpage{089} \fi 
       \ifnum\pagenum=91 \def\dpage{091} \fi \ifnum\pagenum=93 \def\dpage{093} \fi \ifnum\pagenum=95 \def\dpage{095} \fi 
       \ifnum\pagenum=97 \def\dpage{097} \fi \ifnum\pagenum=99 \def\dpage{099} \fi}
\def\reference{ \renewcommand{\baselinestretch}{1} \normalsize
                \renewcommand\section{\@startsection {section}{1}{0mm}%
             {-1.2\baselineskip}%
             {0.5\baselineskip}%
             {\raggedright \@kesect }}
          \section*{\refname}
     \list{}
       {\setlength\leftmargin{2em}%{2.3em}
        \setlength\labelsep{0.3em}
        \setlength\itemindent{-2em}
        \parsep  0pt\relax
        \itemsep\parsep
        \usecounter{enumiv}%
        \let\p@enumiv\@empty
        \def\theenumiv{\arabic{enumiv}}%
        }%
  \let\newblock\bibitem\relax%
}
\def\endreference{%\endenumerate
      \endlist
      \ifx\@submit\@final %
      \samepage
      \begin{center}
      \par \smallskip
           \if@korean
               \footnotesize{\sff \@receiveyear 년 \@receivemonth 월  \@receiveday 일 접수; \@reviseyear 년 \@revisemonth 월 \@reviseday 일 수정; \@acceptyear 년 \@acceptmonth 월 \@acceptday 일 채택}
               \label{sec:\theendpage}
            \else
               \footnotesize{\textsl{Received \convert{\@receivemonth} \@receiveday, \@receiveyear\/; Revised \convert{\@revisemonth} \@reviseday,  \@reviseyear\/; Accepted \convert{\@acceptmonth} \@acceptday, \@acceptyear}}
                  \label{sec:\theendpage}
            \fi
      \end{center}
      \else \fi
}

\def\footnoterule{\kern-3\p@
  \hrule \@width 6cm \kern 2.6\p@} % the \hrule is .4pt high

%       ****************************************
%       *           ADDRESS MACROS             *
%       ****************************************
%
\def\@makeenmark{$^{\@theenmark}$}

\newdimen\endnotesep
\newcommand{\same}[1][]{\hskip 5pt \leavevmode\llap{\hbox{$^{#1}$}\nolinebreak[4]}}%
\newcommand{\address}[1][]{\xdef\@theenmark{#1}\@endnotemark\@endnotetext}
\let\@doanenote=0
\let\@doanenotes=0
\let\@endanenote=0

\newwrite\@enotes
\newif\if@enotesopen \global\@enotesopenfalse

\def\@openenotes{\immediate\openout\@enotes=\jobname.ent\relax
      \global\@enotesopentrue}

\long\def\@endnotetext#1{%
     \if@enotesopen \else \@openenotes \fi
     \immediate\write\@enotes{\@doanenote{\@theenmark}}%
     \begingroup
        \def\next{\mbox{#1{ }}}%
        \newlinechar='40
        \immediate\write\@enotes{\meaning\next}%
     \endgroup
     \immediate\write\@enotes{\@endanenote}}

\def\endnotemark{\xdef\@theenmark{\theendnote}\@endnotemark}

\def\@endnotemark{\leavevmode\ifhmode
  \edef\@x@sf{\the\spacefactor}\fi\@makeenmark
   \ifhmode\spacefactor\@x@sf\fi\relax}

\def\enoteformat{%\rightskip\z@ \leftskip\z@ \parindent=0em%1.8em
     \leavevmode\llap{\hbox{$^{\@theenmark}$}\nolinebreak[4]}}
\def\enotesize{\footnotesize}

\def\theendnotes{\immediate\closeout\@enotes \global\@enotesopenfalse
  \begingroup
    \makeatletter
    \def\@doanenote##1##2>{\def\@theenmark{##1}%\par\begingroup
        \edef\@currentlabel{\csname p@endnote\endcsname\@theenmark} %DW
        \enoteformat}
    \def\@endanenote{}%\par\endgroup}%
    \def\ETC.{\errmessage{Some long endnotes will be truncated; %
                            use BIGLATEX to avoid this}%
          \def\ETC.{\relax}}
    \enotesize
    \begin{center}
    \vspace*{5pt}
    {\normalsize\sffamily%\fontsize{9}{9}\selectfont
    \@input{\jobname.ent}}%
    \end{center}
  \endgroup}
\endinput
           
           

