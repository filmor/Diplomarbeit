%%% ====================================================================
%%% @LaTeX-file{
%%%    filename  = "ml-gen.cls",
%%% }
%%% ====================================================================
%%% derived from gen-j-l.cls ams
% new options
% arxiv, mldraft, final
% option a4paper can be given
% ams classes then adjust paperwidth, paperheight
% and textheight to 54.5pc

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ml-gen}
             [2011/07/20]

\newcommand{\ml@class}{amsart}
\AtEndOfClass{\let\article@logo\ml@logo}
\DeclareOption{book}{\def\ml@class{amsbook}}
\DeclareOption{proc}{\def\ml@class{amsproc}\AtEndOfClass{\let\procart@logo\ml@logo}}

\def\mlSetRCSInfo#1{\def\mlRCSInfo{#1}}              
\InputIfFileExists{hginfo.blg}{}{\mlSetRCSInfo\relax}
\def\mlSetVersion#1{\def\mlVersion{#1}}%
\let\mlVersion\mlRCSInfo

\def\ml@logo{%
  \set@logo{%
    \publname\ \@dateposted\\
    Version: \mlVersion
    }}    

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\ml@class}}

\DeclareOption{arxiv}{%
\AtEndOfClass{
\def\publname{Posted on ArXiV}
\def\datename{\textit{Version}:}
\textwidth33pc
\calclayout}

\PassOptionsToClass{twoside,11pt}{\ml@class}
\PassOptionsToPackage{final}{math60}
}

\DeclareOption{preprint}{%
\AtEndOfClass{
\def\publname{PREPRINT}
\def\datename{\textit{Version}:}
\textwidth33pc
\calclayout}

\PassOptionsToClass{twoside,11pt}{\ml@class}
\PassOptionsToPackage{final}{math60}
}

\DeclareOption{draft}{%
\PassOptionsToClass{oneside,draft,11pt}{\ml@class}
\AtEndOfClass{%
    \def\publname{DRAFT}
    \def\datename{\textit{Versioninfo}:}
    \textwidth33pc
    \calclayout
    \advance \oddsidemargin-1in
    \advance \evensidemargin-1in
    \advance \topmargin-1in
    \advance \marginparwidth2cm
    }
}

\DeclareOption{submit}{%
\PassOptionsToClass{twoside,final}{\ml@class}
    \AtEndOfClass{
        \def\datename{\textit{Date}:}}
\PassOptionsToPackage{final}{math60}
}

\ProcessOptions
\LoadClass{\ml@class}[1996/10/24]
\RequirePackage{math60}             

% copied from amsclass.cls
% only title font changed
% consider latex packages titlesec, sectsty, fncychap, memoir ??
\def\part{\@startsection{part}{0}%
\z@{\linespacing\@plus\linespacing}{.5\linespacing}%
{\Large\bfseries\sffamily\raggedright}}

\def\section{\@startsection{section}{1}%
  \z@{0.7\linespacing\@plus\linespacing}{.5\linespacing}%
  {\large\bfseries\sffamily}}

\def\subsection{\@startsection{subsection}{2}%
  \z@{.5\linespacing\@plus.7\linespacing}{-.5em}%
  {\normalfont\bfseries\sffamily}}
  
\def\subsubsection{\@startsection{subsubsection}{3}%
 \z@{.5\linespacing\@plus.7\linespacing}{-.5em}%
{\normalfont\sffamily\itshape}}
 
%\def\@secnumfont{\bfseries}  
%\def{\datename}{\textit{Version}:}
\def\ISSN{\relax}

\endinput
