\name{brain_extraction}
\alias{brain_extraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brain Extraction
}
\description{
Applies a brain mask to the BOLD image data.
}
\usage{
brain_extraction(
  bold_img,
  mask_img,
  return_nifti = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bold_img}{
NIfTI object for the BOLD data.
}
  \item{mask_img}{
NIfTI object for the brain mask.
}
  \item{return_nifti}{
Return a NIfTI object if TURE, otherwise return a numerical array. (default = FALSE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NIfTI object containing the masked BOLD data.
}
\references{
[1] img_data
https://www.rdocumentation.org/packages/oro.nifti/versions/0.11.4/topics/img_data-methods

[2] as.nifti
https://www.rdocumentation.org/packages/oro.nifti/versions/0.11.4/topics/as.nifti
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
bold_brain <- brain_extraction(bold_img, mask_img)
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
