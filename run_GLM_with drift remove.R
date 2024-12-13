# run script for testing the GLM and modified GLM functions
# loads a minimal fmri dataset and runs the GLM and modified GLM functions
# task used is a motor task, with 2mm and 3mm resolution data


# Install and load packages
if (!requireNamespace("oro.nifti", quietly = TRUE)) {
  install.packages("oro.nifti")
}
if (!requireNamespace("signal", quietly = TRUE)) {
  install.packages("signal")
}
if(!require('caret')) {
  install.packages('caret')
}


library(oro.nifti)
library(signal)
library(caret)

base_dir <- "/Users/hanyang/Desktop/Umich/BIOSTAT 615/Project"

## Step 1: Load the data, perform brain extraction using pre-calculated masks
# Load the local NIfTI images
# t1w_file <- file.path(base_dir,"sub-DF_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii")
bold2mm_file <- file.path(base_dir,"sub-DF_task-motor2mm_space-MNI152NLin2009cAsym_desc-preproc_bold.nii")
# bold3mm_file <- file.path(base_dir,"sub-DF_task-motor3mm_space-MNI152NLin2009cAsym_desc-preproc_bold.nii")
mask2mm_file <- file.path(base_dir,"sub-DF_task-motor2mm_space-MNI152NLin2009cAsym_desc-brain_mask.nii")
# mask3mm_file <- file.path(base_dir,"sub-DF_task-motor3mm_space-MNI152NLin2009cAsym_desc-brain_mask.nii")

# t1w_img <- readNIfTI(t1w_file, reorient = FALSE) # structural image
bold2mm_img <- readNIfTI(bold2mm_file, reorient = FALSE) # bold 2mm resolution
# bold3mm_img <- readNIfTI(bold3mm_file, reorient = FALSE) # bold 3mm resolution
mask2mm_img <- readNIfTI(mask2mm_file, reorient = FALSE) # brain mask 2mm resolution
# mask3mm_img <- readNIfTI(mask3mm_file, reorient = FALSE) # brain mask 3mm resolution

# extract the bold data from the 4D image
bold2mm <- img_data(bold2mm_img)
# bold3mm <- img_data(bold3mm_img)
mask2mm <- img_data(mask2mm_img)
# mask3mm <- img_data(mask3mm_img)

# brain extraction of fmri using brain masks
mask2mm_4d <- array(mask2mm, dim = dim(bold2mm_img))  # Expand the mask to 4D
# mask3mm_4d <- array(mask3mm, dim = dim(bold3mm_img)) 

# Apply the brain mask by element-wise multiplication
bold2mm_brain <- bold2mm_img * mask2mm_4d
# bold3mm_brain <- bold3mm_img * mask3mm_4d

# reshape the bold data to a 2D matrix
# bold2mm_data_2d <- matrix(bold2mm_data, nrow = n_vols, byrow = TRUE)
# bold3mm_data_2d <- matrix(bold3mm_data, nrow = n_vols, byrow = TRUE)

# orthographic(t1w_img, main = "Orthographic View of T1w Image")

## Step 2: Convolve task timing with HRF to create task regressor b

# Load data
task_file_1 <- file.path(base_dir, "task-motor_events.tsv")
task_data <- read.table(task_file_1, sep = "\t", header = TRUE)
task_file_2 <- file.path(base_dir, "sub-DF_task-motor2mm_desc-confounds_timeseries.tsv")
confounds <- read.table(task_file_2, sep = "\t", header = TRUE)
head(task_data)
head(confounds)

# Processing the confounds data

# Impute NAs with column means
confounds[confounds == "n/a"] <- NA

# Convert all columns to numeric
confounds <- as.data.frame(lapply(confounds, function(col) {
  as.numeric(as.character(col))
}))

confounds_imputed <- confounds  # Create a copy of the data frame

for (col in names(confounds_imputed)) {
  confounds_imputed[[col]][is.na(confounds_imputed[[col]])] <- 
    mean(as.numeric(confounds_imputed[[col]]), na.rm = TRUE)
}

# Verify that there are no more NAs
cat("Number of remaining NAs:", sum(is.na(confounds_imputed)), "\n")

# Remove columns with low variance
confounds_df_filtered <- confounds_imputed[, sapply(confounds_imputed, function(col) sd(col, na.rm = TRUE) > 1e-6)]

# Compute the correlation matrix
cor_matrix <- cor(confounds_df_filtered, use = "pairwise.complete.obs")

# Find highly correlated columns
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8)

# Drop highly correlated columns
confounds_df_selected <- confounds_df_filtered[, -high_corr]

# Get names of selected columns
selected_columns <- colnames(confounds_df_selected)

# Print the number of selected columns and their names
cat("Number of selected columns:", length(selected_columns), "\n")
print(selected_columns)

# Save the filtered confounds to a new data frame
design_matrix_confounds <- confounds_df_selected

TR <- bold2mm_img@pixdim[4]  # Repetition time in seconds
n_vols <- dim(bold2mm_img)[4]  # Number of volumes in the 4D fMRI data
task_regressor <- rep(0, n_vols) # Initialize task regressor vector

# Assume task_data has columns 'onset' and 'duration' in seconds
for (i in 1:nrow(task_data)) {
  onset <- round(task_data$onset[i] / TR)
  duration <- round(task_data$duration[i] / TR)
  task_regressor[onset:(onset + duration - 1)] <- 1
}

# Perform convolutions
spm_hrf <- function(TR, p = c(6, 16, 1, 1, 6, 0, 32)) {
  dt <- TR / 16
  u <- seq(0, p[7] / dt) * dt
  hrf <- dgamma(u, p[1] / p[3], dt / p[3]) - dgamma(u, p[2] / p[4], dt / p[4]) / p[5]
  hrf <- hrf / max(hrf) # normalize to max of 1
  return(hrf)
}

hrf <- spm_hrf(TR)
task_conv <- convolve(task_regressor, rev(hrf), type = "open")  # Perform convolution
task_conv_trimmed <- tail(task_conv, 114)
plot(task_conv, type = "l", main = "HRF-convolved Task Regressor", xlab = "Time (TR)", ylab = "Amplitude") # Check regressor

# Combine task_conv and confounds to get the design matrix
task_regressor_matrix <- matrix(task_conv_trimmed, ncol = 1)  # Convert to a matrix

# Combine task regressors with selected confounds
design_matrix <- cbind(Intercept = 1, task_regressor_matrix, design_matrix_confounds)

# Check the dimensions of the final design matrix
cat("Design matrix dimensions:", dim(design_matrix), "\n")

## step 3: manually create the contrast vector based on trial types
contrast_vec <- array(c(-1, 1), 141)

# task_conv is 'b' in the GLM equation: yi=βbi+fi+ni,i=1,…,N
# where yi is the observed fMRI signal at time i, fi is the drift, and ni is the noise

# Perform GLM analysis
glm_voxelwise <- function(Y, design_matrix, contrast_vec, lambda = 1e4, block_size = 10) {
  # Step 1: Convert design_matrix (list) to matrix
  X <- do.call(cbind, design_matrix)  # Combine list elements into a matrix
  contrast_vector <- as.numeric(contrast_vec)
  
  # Validate dimensions
  if (length(contrast_vector) != ncol(X)) {
    stop("Contrast vector length must match the number of columns in the design matrix")
  }
  
  # Dimensions of input data
  x_dim <- dim(Y)[1]
  y_dim <- dim(Y)[2]
  z_dim <- dim(Y)[3]
  n <- dim(Y)[4]
  k <- ncol(X)
  
  # Precompute the regularized pseudo-inverse of X
  XtX <- t(X) %*% X
  XtX_inv <- tryCatch({
    solve(XtX)  # Standard GLM
  }, error = function(e) {
    solve(XtX + lambda * diag(k))  # Ridge regression
  })
  XtX_inv_Xt <- XtX_inv %*% t(X)
  
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
    
    # Solve for β coefficients for valid voxels
    beta_block <- matrix(0, nrow = prod(dim(Y_block)[1:3]), ncol = k)
    beta_block[valid_voxels, ] <- XtX_inv_Xt %*% t(Y_block_flat[valid_voxels, ])
    
    # Compute contrast image for the block
    contrast_block <- beta_block %*% contrast_vector
    
    # Reshape back to the original dimensions
    beta_block_reshaped <- array(beta_block, dim = c(dim(Y_block)[1], dim(Y_block)[2], dim(Y_block)[3], k))
    contrast_block_reshaped <- array(contrast_block, dim = c(dim(Y_block)[1], dim(Y_block)[2], dim(Y_block)[3]))
    
    # Store results
    beta_images[, , z:z_end, ] <- beta_block_reshaped
    contrast_image[, , z:z_end] <- contrast_block_reshaped
  }
  
  # Return results
  list(beta_images = beta_images, contrast_image = contrast_image)
}


# Usage
beta_results <- glm_voxelwise(bold2mm_brain, design_matrix, contrast_vec)
writeNIfTI(as.nifti(beta_results$contrast_image), "contrast_image")

# Load the NIfTI file
nii_file <- readNIfTI("contrast_image.nii")

# View orthogonal slices
orthographic(nii_file, main = "Contrast Image")


#Donglin Liu
## Step 4: Implement Modified GLM in Wavelet Domain
library(wavelets)

# Function to perform Modified GLM in Wavelet Domain
modified_glm_wavelet <- function(Y, design_matrix, drift_levels = 3) {
  # Step 1: Apply Discrete Wavelet Transform (DWT) to the data
  dwt_results <- dwt(Y, filter = "haar", boundary = "periodic")  # Using Haar wavelet
  
  # Extract coefficients at different levels
  wavelet_coefficients <- dwt_results@W
  scaling_coefficients <- dwt_results@V
  
  # Step 2: Remove drift by zeroing coefficients at coarse levels
  for (level in seq_along(wavelet_coefficients)) {
    if (level > drift_levels) {
      wavelet_coefficients[[level]] <- 0  # Ignore finer levels for drift
    }
  }
  
  # Step 3: Perform Bayesian estimation in wavelet domain
  # Using ridge regression-like Bayesian estimation
  X <- design_matrix
  XtX <- t(X) %*% X
  lambda <- 1e4  # Regularization parameter
  XtX_inv <- solve(XtX + lambda * diag(ncol(X)))
  XtX_inv_Xt <- XtX_inv %*% t(X)
  
  beta_estimates <- lapply(wavelet_coefficients, function(coeff) {
    if (!is.null(coeff)) {
      coeff_flat <- as.vector(coeff)
      valid_voxels <- !is.na(coeff_flat)
      beta <- rep(NA, length(coeff_flat))
      beta[valid_voxels] <- XtX_inv_Xt %*% coeff_flat[valid_voxels]
      matrix(beta, nrow = dim(coeff))
    } else {
      NULL
    }
  })
  
  # Step 4: Reconstruct the signal without drift
  dwt_results@W <- wavelet_coefficients  # Use modified coefficients
  reconstructed_signal <- idwt(dwt_results)  # Perform Inverse DWT
  
  # Return the reconstructed signal and estimated beta values
  list(reconstructed_signal = reconstructed_signal, beta_estimates = beta_estimates)
}

# Apply Modified GLM in Wavelet Domain
modified_results <- modified_glm_wavelet(bold2mm_brain, design_matrix, drift_levels = 3)

# Save reconstructed signal as NIfTI file
reconstructed_signal <- modified_results$reconstructed_signal
writeNIfTI(as.nifti(reconstructed_signal), "reconstructed_signal")

# Load and visualize reconstructed signal
nii_reconstructed <- readNIfTI("reconstructed_signal.nii")
orthographic(nii_reconstructed, main = "Reconstructed Signal")

# Save beta coefficients for visualization
beta_wavelet <- modified_results$beta_estimates
saveRDS(beta_wavelet, "beta_wavelet_coefficients.rds")

# Plot Beta Coefficients for Key Wavelet Levels
for (level in seq_along(beta_wavelet)) {
  if (!is.null(beta_wavelet[[level]])) {
    plot(beta_wavelet[[level]], type = "l", main = paste("Beta Coefficients at Wavelet Level", level),
         xlab = "Voxel Index", ylab = "Beta Value")
  }
}
