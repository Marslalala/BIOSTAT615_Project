\name{generate_task_regressor}
\alias{generate_task_regressor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Task Regressor (necessary for the design matrix)
}
\description{
Convolves task timing data with the HRF to create a task regressor.
}
\usage{
generate_task_regressor(
  task_file,
  TR
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{task_file}{
Path to the task timing data (TSV format).
}
  \item{TR}{
Repetition time (in seconds) from the BOLD data.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector representing the task regressor.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yang Han

William Tackett

Donglin Liu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
task_regressor <- generate_task_regressor("path/to/task.tsv", TR = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
