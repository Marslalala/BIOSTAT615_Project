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

library(oro.nifti)
library(signal)

base_dir <- "/Users/will/Documents/Documents - Mac (2) 2/UM/Fall24/BIOSTAT615/project/"

## Step 1: Load the data, perform brain extraction using pre-calculated masks
# Load the local NIfTI images
t1w_file <- file.path(base_dir,"data/sub-DF_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz")
bold2mm_file <- file.path(base_dir,"data/sub-DF_task-motor2mm_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
bold3mm_file <- file.path(base_dir,"data/sub-DF_task-motor3mm_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
mask2mm_file <- file.path(base_dir,"data/sub-DF_task-motor2mm_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz")
mask3mm_file <- file.path(base_dir,"data/sub-DF_task-motor3mm_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz")

t1w_img <- readNIfTI(t1w_file, reorient = FALSE) # structural image
bold2mm_img <- readNIfTI(bold2mm_file, reorient = FALSE) # bold 2mm resolution
bold3mm_img <- readNIfTI(bold3mm_file, reorient = FALSE) # bold 3mm resolution
mask2mm_img <- readNIfTI(mask2mm_file, reorient = FALSE) # brain mask 2mm resolution
mask3mm_img <- readNIfTI(mask3mm_file, reorient = FALSE) # brain mask 3mm resolution

# extract the bold data from the 4D image
bold2mm <- img_data(bold2mm_img)
bold3mm <- img_data(bold3mm_img)
mask2mm <- img_data(mask2mm_img)
mask3mm <- img_data(mask3mm_img)

# brain extraction of fmri using brain masks
mask2mm_4d <- array(mask2mm, dim = dim(bold2mm_img))  # Expand the mask to 4D
mask3mm_4d <- array(mask3mm, dim = dim(bold3mm_img)) 

# Apply the brain mask by element-wise multiplication
bold2mm_brain <- bold2mm_img * mask2mm_4d
bold3mm_brain <- bold3mm_img * mask3mm_4d

# reshape the bold data to a 2D matrix
# bold2mm_data_2d <- matrix(bold2mm_data, nrow = n_vols, byrow = TRUE)
# bold3mm_data_2d <- matrix(bold3mm_data, nrow = n_vols, byrow = TRUE)

orthographic(t1w_img, main = "Orthographic View of T1w Image")

## Step 2: Convolve task timing with HRF to create task regressor b

# Load data
task_file <- file.path(base_dir, "data/task-motor_events.tsv")
task_data <- read.table(task_file, sep = "\t", header = TRUE)
head(task_data)

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
hrf <- spm_hrf(TR)
task_conv <- convolve(task_regressor, rev(hrf), type = "open")  # Perform convolution
plot(task_conv, type = "l", main = "HRF-convolved Task Regressor", xlab = "Time (TR)", ylab = "Amplitude") # Check regressor

# task_conv is 'b' in the GLM equation: yi=βbi+fi+ni,i=1,…,N
# where yi is the observed fMRI signal at time i, fi is the drift, and ni is the noise

# Perform GLM analysis



