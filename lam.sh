#!/bin/bash
#SBATCH -n 64
#SBATCH --partition=rome
#SBATCH -t 120:00:00

# Script options
cleanup_scratch=false
backup_google=false
fetch_restarts=false

export MODEL_PARAM_FILE=output.dat
export MODEL_NAME="t1_64"
export INPUT_FOLDER=$PWD
export OUTPUT_FOLDER=$INPUT_FOLDER/ModelOutputFiles/$MODEL_NAME
export RESTART_FOLDER=$INPUT_FOLDER/ModelRestartFiles

module load 2023 PETSc/3.20.3-foss-2023a HDF5/1.14.0-gompi-2023a h5py/3.9.0-foss-2023a SWIG/4.1.1-GCCcore-12.3.0 Ninja/1.11.1-GCCcore-12.3.0 

# Create directories if they don't exist
mkdir -p "$OUTPUT_FOLDER"
mkdir -p "$RESTART_FOLDER"

# Set error and output file names
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
error_file="${MODEL_NAME}_${timestamp}.err"
output_file="${MODEL_NAME}_${timestamp}.out"

# Redirect stdout and stderr to output and error files
exec > >(tee "$output_file") 2> >(tee "$error_file" >&2)

function print_script_options() {
    echo "UW_VERSION: $IMAGE_VERSION"
    echo "Script options:"
    echo "MODEL_NAME: $MODEL_NAME"
    echo "cleanup_scratch: $cleanup_scratch"
    echo "backup_google: $backup_google"
    echo "fetch_restarts: $fetch_restarts"
    echo "MODEL_SCRIPT: $MODEL_SCRIPT"
    echo "INPUT_FOLDER: $INPUT_FOLDER"
    echo "OUTPUT_FOLDER: $OUTPUT_FOLDER"
    echo "RESTART_FOLDER: $RESTART_FOLDER"
    echo "----------------------------------"
}

function fetch_restart_files() {
    if $fetch_restarts; then
        echo "Fetching restart files"
        . fetchRestarts.sh >"$MODEL_NAME"_restartCopy_log.out
        wait
    fi
}

function backup_to_google_drive() {
    if $backup_google; then
        cuDir="$PWD"
        echo "Backing up to Google Drive"
        ModTMP="${MODEL_NAME}_TEMP_FOLDER"
        ls $(cat $ModTMP)
        cd $(cat $ModTMP)
        echo $PWD $cuDir
        rclone copy ./ gdalk:Snellius3D/SAnRef/ -vvv
        cd $cuDir
        wait
    fi
}

function create_temp_directory() {
    SCRATCH_SHARED_DIR=$(mktemp -d -p /scratch-shared)
    echo "${SCRATCH_SHARED_DIR}" >"$MODEL_NAME"_TEMP_FOLDER
    echo "Created temporary directory: $SCRATCH_SHARED_DIR"
}

function copy_scripts_and_checkpoint_files() {
    echo "Copying scripts and checkpoint files to temporary directory..."
    cp -v "$INPUT_FOLDER"/* "$SCRATCH_SHARED_DIR"
    cp -v -r "$INPUT_FOLDER"/markers "$SCRATCH_SHARED_DIR/markers"

    cp -v -r "$RESTART_FOLDER/$MODEL_NAME" "$SCRATCH_SHARED_DIR"
    echo "Finished copying scripts and checkpoint files to temporary directory."
}

function list_temp_directory() {
    echo "Listing temporary directory:"
    cd "$SCRATCH_SHARED_DIR"
    ls -l -r
    echo "Listing temporary directory Model Files:"
    ls -l -r $SCRATCH_SHARED_DIR/$MODEL_NAME
    echo "----------------------------------"
}

function load_modules_and_run_model() {
    echo "********************************** Model Running ************************************"
    srun -n 64 /home/alaik/lam/LaMEM/bin/opt/LaMEM -ParamFile $MODEL_PARAM_FILE
}

function copy_output_files() {
    rsync -a -v --ignore-existing --exclude '*.sif' "$PWD/" "$OUTPUT_FOLDER/"
    rsync -a -v --include '*.log' "$PWD/" "$OUTPUT_FOLDER/"
    echo "Output copied to $OUTPUT_FOLDER"
}

function cleanup_temp_directory() {
    if $cleanup_scratch; then
        rm -rf "$SCRATCH_SHARED_DIR"
        echo "Temporary directory $SCRATCH_SHARED_DIR removed"
    fi
}

# Main script
print_script_options
fetch_restart_files
backup_to_google_drive
create_temp_directory
copy_scripts_and_checkpoint_files
list_temp_directory
load_modules_and_run_model
copy_output_files
cleanup_temp_directory
