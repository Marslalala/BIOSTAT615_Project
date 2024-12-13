\name{glm_voxelwise}
\alias{glm_voxelwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Linear Model (GLM) for voxel-wise Analysis
}
\description{
Fits a GLM to fMRI data at each voxel and computes a contrast image.
}
\usage{
glm_voxelwise(
  Y,
  design_matrix,
  contrast_vector,
  lambda = 1e4,
  block_size = 10
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
4D array of fMRI data (masked BOLD signal).
}
  \item{design_matrix}{
Matrix representing the design (task regressor + confounds, if any).
}
  \item{contrast_vector}{
Numeric vector specifying the contrast of interest.
}
  \item{lambda}{
Regularization parameter for ridge regression (default = 1e4).
}
  \item{block_size}{
Number of slices to process simultaneously (default = 10).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing:
  \describe{
     \item{beta_images}{4D array of estimated beta coefficients.}
     \item{contrast_image}{3D array of the contrast image.}
}
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
results <- glm_voxelwise(bold_brain, task_regressor_matrix, contrast_vector)
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
