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
glm_voxelwise <- function(Y, design_mat, contrast_vec) {
  # # Step 1: Load the design matrix and contrast vector
  # X <- as.matrix(read.table(design_mat))
  # contrast_vector <- as.numeric(scan(contrast_vec, what = double()))
  # 
  # # Step 2: Load images into a 4D array
  # # Read all images and stack them into a 4D array
  # image_list <- lapply(image_files, readNIfTI)
  # Y <- array(unlist(image_list), dim = c(dim(image_list[[1]]), length(image_list)))
  
  # Step 1: Convert design_matrix (list) to matrix and validate inputs
  X <- do.call(cbind, design_matrix)  # Combine list elements into a matrix
  contrast_vector <- as.numeric(contrast_vec)  # Ensure contrast_vec is numeric
  
  # Validate dimensions
  if (length(contrast_vector) != ncol(X)) {
    stop("Contrast vector length must match the number of columns in the design matrix")
  }
  
  # Step 3: Dimensions of input data
  x_dim <- dim(Y)[1]
  y_dim <- dim(Y)[2]
  z_dim <- dim(Y)[3]
  n <- dim(Y)[4]
  k <- ncol(X) # number of columns in X (number of factors)
  lambda <- 1e4
  
  # Step 4: Create an empty 4D array for β coefficient images
  beta_images <- array(0, dim = c(x_dim, y_dim, z_dim, k))
  
  # Step 5: Loop over each voxel position (i, j, l)
  for (i in 1:x_dim) {
    for (j in 1:y_dim) {
      for (l in 1:z_dim) {
        
        # Step 5a: Extract voxel values across images for current voxel position
        voxel_values <- Y[i, j, l, ]
        
        # Step 5b: Solve for β coefficients at this voxel position
        if (all(!is.na(voxel_values) & !is.nan(voxel_values))) { # Check for valid values
          # Attempt standard GLM solution
          beta <- tryCatch({
            solve(t(X) %*% X) %*% t(X) %*% voxel_values  # Standard GLM
          }, error = function(e) {
            # If singular, apply ridge regression
            solve(t(X) %*% X + lambda * diag(k)) %*% t(X) %*% voxel_values
          })
          
          # Store each β coefficient in the created beta_images array
          for (factor in 1:k) {
            beta_images[i, j, l, factor] <- beta[factor]
          }
        }
      }
    }
  }
  
  # Step 6: Apply the contrast vector to compute the contrast image
  contrast_image <- array(0, dim = c(x_dim, y_dim, z_dim))
  for (factor in 1:k) {
    contrast_image <- contrast_image + contrast_vector[factor] * beta_images[, , , factor]
  }
  
  # Step 7: Return the list of β coefficient images and the contrast image
  list(beta_images = beta_images, contrast_image = contrast_image)
}

# Example usage:
beta_results <- glm_voxelwise(bold2mm_brain, design_matrix, contrast_vec)
writeNIfTI(beta_results$contrast_image, "contrast_image.nii")

