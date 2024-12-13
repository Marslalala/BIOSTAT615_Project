\name{save_glm_results}
\alias{save_glm_results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Save GLM Results to NIfTI Files
}
\description{
Saves contrast image and beta coefficient maps from GLM results as NIfTI files.
}
\usage{
save_glm_results(
  contrast_image,
  beta_images,
  output_dir = NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrast_image}{
3D array of contrast values.
}
  \item{beta_images}{
4D array of beta coefficients.
}
  \item{output_dir}{
(Optional) Directory where the results will be saved. If not provided, saves to the current working directory.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
None. Saves NIfTI files to the specified directory or current working directory.
}
\references{
[1] as.nifti
https://www.rdocumentation.org/packages/oro.nifti/versions/0.11.4/topics/as.nifti

[2] writeNIfTI
https://www.rdocumentation.org/packages/oro.nifti/versions/0.11.4/topics/writeNIfTI-methods
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
save_glm_results(contrast_image, beta_images)
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
