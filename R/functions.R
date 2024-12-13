### This file contains all the functions in the package

# 1. Data Loading

#' Load NIfTI Data
#'
#' Reads a NIfTI file and returns the image data. If the file path is relative,
#' it is resolved relative to the current working directory.
#'
#' @param file_path Path to the NIfTI file. Can be relative to the working directory.
#' @return An object containing the NIfTI image data.
#' @examples
#' bold_img <- load_nifti("bold.nii")  # File in working directory
#' bold_img <- load_nifti("/absolute/path/to/bold.nii")  # Absolute path
load_nifti <- function(file_path) {
  if (!requireNamespace("oro.nifti", quietly = TRUE)) {
    stop("Package 'oro.nifti' is required. Please install it.")
  }
  
  # Resolve file path relative to the current working directory
  resolved_path <- normalizePath(file_path, mustWork = FALSE)
  
  # Check if the file exists
  if (!file.exists(resolved_path)) {
    stop("File not found: ", resolved_path)
  }
  
  # Load and return the NIfTI image
  oro.nifti::readNIfTI(resolved_path, reorient = FALSE)
}


# 2. Observed Data Preprocessing

#' Brain Extraction
#'
#' Applies a brain mask to the BOLD image data.
#'
#' @param bold_img NIfTI object for the BOLD data.
#' @param mask_img NIfTI object for the brain mask.
#' @param return_nifti Return a NIfTI object if TURE, otherwise return a numerical array. (default = FALSE)
#' @return NIfTI object containing the masked BOLD data.
#' @examples
#' bold_brain <- brain_extraction(bold_img, mask_img)
brain_extraction <- function(bold_img, mask_img, return_nifti = FALSE) {
  bold_data <- oro.nifti::img_data(bold_img)
  mask_data <- oro.nifti::img_data(mask_img)
  
  if (!all(dim(bold_data)[1:3] == dim(mask_data))) {
    stop("Dimensions of BOLD and mask images do not match!")
  }
  
  # Expand the mask to 4D
  mask_4d <- array(mask_data, dim = dim(bold_data))
  
  # Apply the mask
  masked_bold <- bold_data * mask_4d
  
  if (return_nifti) {
    return(oro.nifti::as.nifti(masked_bold, template = bold_img))
  }
  else {
    return(masked_bold)
  }
}


# 3. Task Regressor Generating

#' Generate Task Regressor (necessary for the design matrix)
#'
#' Convolves task timing data with the HRF to create a task regressor.
#'
#' @param task_file Path to the task timing data (TSV format).
#' @param TR Repetition time (in seconds) from the BOLD data.
#' @param duration Duration of each condition in seconds (default = 16).
#' @return A vector representing the task regressor.
#' @examples
#' task_regressor <- generate_task_regressor("path/to/task.tsv", 2)
generate_task_regressor <- function(task_file, TR) {
  # Load task timing data
  task_data <- read.table(task_file, sep = "\t", header = TRUE)
  
  # Define the parameters for hrf
  hrf_params = c(6, 16, 1, 1, 6, 0, 32)
  
  # Define SPM-style HRF
  spm_hrf <- function(TR, P) {
    t <- seq(0, P[7] / TR, by = TR)
    hrf <- dgamma(t, shape = P[1] / P[3], scale = P[3]) - 
      dgamma(t, shape = P[2] / P[4], scale = P[4]) / P[5]
    hrf <- hrf / sum(hrf)
    return(hrf)
  }
  
  # Create HRF
  hrf <- spm_hrf(TR, hrf_params)
  
  # Generate binary task regressor from trial types
  n_vols <- max(task_data$onset + task_data$duration) / TR
  n_reps <- task_data$duration[1] / 2
  output_vector <- unlist(lapply(task_data$trial_type, function(tt) {
    if (tt == "stimulus") {
      rep(1, n_reps)
    } else if (tt == "baseline") {
      rep(0, n_reps)
    } else {
      stop("Unknown trial_type")
    }
  }))
  
  # Convolve with HRF
  task_regressor <- convolve(output_vector, rev(hrf), type = "filter")
  
  # Ensure regressor length matches expected volumes
  task_regressor <- task_regressor[1:n_vols]
  return(task_regressor)
}


# 4. GLM Function

#' General Linear Model (GLM) for voxel-wise Analysis
#'
#' Fits a GLM to fMRI data at each voxel and computes a contrast image.
#'
#' @param Y 4D array of fMRI data (masked BOLD signal).
#' @param design_matrix Matrix representing the design (task regressor + confounds, if any).
#' @param contrast_vector Numeric vector specifying the contrast of interest.
#' @param lambda Regularization parameter for ridge regression (default = 1e4).
#' @param block_size Number of slices to process simultaneously (default = 10).
#' @return A list containing:
#'   \describe{
#'     \item{beta_images}{4D array of estimated beta coefficients.}
#'     \item{contrast_image}{3D array of the contrast image.}
#'   }
#' @examples
#' results <- glm_voxelwise(bold_brain, task_regressor_matrix, c(-1, 1))
glm_voxelwise <- function(Y, design_matrix, contrast_vector, lambda = 1e4, block_size = 10) {
  # Validate inputs
  if (length(contrast_vector) != ncol(design_matrix)) {
    stop("Contrast vector length must match the number of columns in the design matrix.")
  }
  
  # Dimensions of input data
  x_dim <- dim(Y)[1]
  y_dim <- dim(Y)[2]
  z_dim <- dim(Y)[3]
  n <- dim(Y)[4]
  k <- ncol(design_matrix)
  
  # Precompute the regularized inverse of X
  XtX <- t(design_matrix) %*% design_matrix
  XtX_inv <- tryCatch({
    solve(XtX)  # Standard inverse
  }, error = function(e) {
    solve(XtX + lambda * diag(k))  # Regularized inverse
  })
  XtX_inv_Xt <- XtX_inv %*% t(design_matrix)
  
  # Initialize arrays for results
  beta_images <- array(0, dim = c(x_dim, y_dim, z_dim, k))
  contrast_image <- array(0, dim = c(x_dim, y_dim, z_dim))
  
  # Process data in blocks
  for (z in seq(1, z_dim, by = block_size)) {
    z_end <- min(z + block_size - 1, z_dim)
    
    # Extract the block of data
    Y_block <- Y[, , z:z_end, ]
    Y_block_flat <- matrix(Y_block, nrow = prod(dim(Y_block)[1:3]), ncol = n)
    
    # Skip empty or invalid data
    valid_voxels <- rowSums(!is.na(Y_block_flat)) == n
    
    # Solve for beta coefficients for valid voxels
    beta_block <- matrix(0, nrow = prod(dim(Y_block)[1:3]), ncol = k)
    beta_block[valid_voxels, ] <- XtX_inv_Xt %*% t(Y_block_flat[valid_voxels, ])
    
    # Compute contrast image for the block
    contrast_block <- beta_block %*% contrast_vector
    
    # Reshape back to the original dimensions
    beta_block_reshaped <- array(beta_block, dim = c(dim(Y_block)[1], dim(Y_block)[2], dim(Y_block)[3], k))
    contrast_block_reshaped <- array(contrast_block, dim = c(dim(Y_block)[1], dim(Y_block)[2], dim(Y_block)[3]))
    
    beta_images[, , z:z_end, ] <- beta_block_reshaped
    contrast_image[, , z:z_end] <- contrast_block_reshaped
  }
  
  list(beta_images = beta_images, contrast_image = contrast_image)
}


# 5. Save the results

#' Save GLM Results to NIfTI Files
#'
#' Saves contrast image and beta coefficient maps from GLM results as NIfTI files.
#'
#' @param contrast_image 3D array of contrast values.
#' @param beta_images 4D array of beta coefficients.
#' @param output_dir (Optional) Directory where the results will be saved. If not provided, saves to the current working directory.
#' @param template_img NIfTI object used as a template for metadata.
#' @return None. Saves NIfTI files to the specified directory or current working directory.
#' @examples
#' save_glm_results(contrast_image, beta_images, template_img = bold_img)
save_glm_results <- function(contrast_image, beta_images, output_dir = NULL) {
  # Use current working directory if output_dir is NULL
  if (is.null(output_dir)) {
    output_dir <- getwd()
  } else {
    # Resolve output directory to absolute path
    output_dir <- normalizePath(output_dir, mustWork = FALSE)
  }
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the contrast image
  contrast_nifti <- oro.nifti::as.nifti(contrast_image)
  contrast_file <- file.path(output_dir, "contrast_image")
  oro.nifti::writeNIfTI(contrast_nifti, contrast_file)
  cat("Contrast image saved to:", contrast_file, "\n")
  
  # Save beta coefficient maps
  for (i in seq_len(dim(beta_images)[4])) {
    beta_nifti <- oro.nifti::as.nifti(beta_images[, , , i])
    beta_file <- file.path(output_dir, paste0("beta_image_", i, ".nii"))
    oro.nifti::writeNIfTI(beta_nifti, beta_file)
    cat("Beta image", i, "saved to:", beta_file, "\n")
  }
}