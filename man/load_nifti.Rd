\name{load_nifti}
\alias{load_nifti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load NIfTI Data
}
\description{
Reads a NIfTI file and returns the image data. If the file path is relative, it is resolved relative to the current working directory.
}
\usage{
load_nifti(file_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path}{
Path to the NIfTI file. Can be relative to the working directory.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object containing the NIfTI image data.
}
\references{
readNIfTI
https://www.rdocumentation.org/packages/dcemriS4/versions/0.20.2/topics/readNIfTI
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
bold_img <- load_nifti("bold.nii")  # File in working directory
bold_img <- load_nifti("/absolute/path/to/bold.nii")  # Absolute path
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
