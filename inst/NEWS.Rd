\name{NEWS}
\title{ppcc News}
\encoding{UTF-8}
\section{Changes in version 1.2}{
  \subsection{Others}{
    \itemize{
      \item Calls to \sQuote{R_registerRoutines} and 
      \sQuote{R_useDynamicSymbols} were re-written.
      \item NAMESPACE was modified towards 
      \sQuote{useDynLib('ppcc', .registration = TRUE, .fixes = "C_")}
    }
  }
}
\section{Changes in version 1.1}{
  \subsection{Bugfixes}{
    \itemize{
      \item The quantile calculation for Pearson type 3
      distribution if shape > 0 was corrected.
      }
   }
}
