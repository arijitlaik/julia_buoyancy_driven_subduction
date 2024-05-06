#!/bin/bash

# Default values
LOGFILE="./pushtohpc.log"
WORKSPACE_DIR="./"  
DESTINATION="alaik@snellius.surf.nl:lam"
DEFAULT_SCRIPT="test.jl"
IS_DRY_RUN=0
IS_SKIP_JULIA=0
IS_SKIP_RSYNC=0
SCRIPT_TO_RUN="$DEFAULT_SCRIPT"

# Function to print usage
print_usage() {
    echo "Usage: $0 [-h] [-n] [-j] [-r] [-s script] [-d directory] [-t target]"
    echo "  -h            Display help"
    echo "  -n            Dry run: display what would be done"
    echo "  -j            Skip running Julia script"
    echo "  -r            Skip performing rsync"
    echo "  -s script     Julia script to run (default: $DEFAULT_SCRIPT)"
    echo "  -d directory  Directory to rsync (default: $WORKSPACE_DIR)"
    echo "  -t target     Target for rsync (default: $DESTINATION)"
    return 1
}

# Function to log messages
log_message() {
    local message="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" >> $LOGFILE
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" 

}

# Function to handle errors
handle_error() {
    local exit_code="$1"
    local action="$2"
    if [ $exit_code -ne 0 ]; then
        log_message "Error during $action. Exiting."
        exit 1
    fi
    log_message "$action completed successfully."
}

# Function to run Julia script
run_julia_script() {
    local script="$1"
    if [ $IS_DRY_RUN -eq 1 ]; then
        log_message "Dry run: Would run Julia script: $script"
    else
        log_message "Running Julia script: $script"
        julia "$script"
        handle_error $? "running Julia script"
    fi
}

# Function to perform rsync
perform_rsync_operation() {
    local dir="$1"
    local dest="$2"
    if [ $IS_DRY_RUN -eq 1 ]; then
        log_message "Dry run: Would perform rsync operation from $dir to $dest"
        rsync -avz --info=progress2 --include 'm*/' --include '*.dat' --include 'Proc*.bin' --include '*.vts' --include '*.sh' --exclude '*' --exclude '.*' --update --dry-run "$dir" "$dest"

    else
        log_message "Starting rsync..."
        # rsync -avz --info=progress2 --include 'm*/'  --include '*.dat'  --include 'Proc*.bin' --include '*.vts' --include '*.sh' --exclude '*' --exclude '.*' "$dir" "$dest"
        rsync -avz --info=progress2 --include 'm*/' --include '*.dat' --include 'Proc*.bin' --include '*.vts' --include '*.sh' --exclude '*' --exclude '.*' --update "$dir" "$dest"

        handle_error $? "rsync operation"
    fi
}

# Function to ssh and run script
ssh_and_run_script() {
    if [ $IS_DRY_RUN -eq 1 ]; then
        log_message "Dry run: Would SSH into server and run script"
    else
        log_message "SSH into the server and running script..."
        ssh -t alaik@snellius.surf.nl "cd lam && sbatch lamsa.sh"
        handle_error $? "SSH or running script"
    fi
}

# Function to parse command line options
parse_command_line_options() {
    local _is_dry_run=$IS_DRY_RUN
    local _is_skip_julia=$IS_SKIP_JULIA
    local _is_skip_rsync=$IS_SKIP_RSYNC
    local _script_to_run=$SCRIPT_TO_RUN
    local _workspace_dir=$WORKSPACE_DIR
    local _destination=$DESTINATION

    while getopts "hnjrs:d:t:" opt; do
        case ${opt} in
            h)
                print_usage
                return 1
                ;;
            n)
                _is_dry_run=1
                ;;
            j)
                _is_skip_julia=1
                ;;
            r)
                _is_skip_rsync=1
                ;;
            s)
                _script_to_run="$OPTARG"
                ;;
            d)
                _workspace_dir="$OPTARG"
                ;;
            t)
                _destination="$OPTARG"
                ;;
            \?)
                echo "Invalid option: -$OPTARG" 1>&2
                print_usage
                return 1
                ;;
            :)
                echo "Option -$OPTARG requires an argument." 1>&2
                print_usage
                return 1
                ;;
        esac
    done
    shift $((OPTIND -1))

    echo "$_is_dry_run $_is_skip_julia $_is_skip_rsync $_script_to_run $_workspace_dir $_destination"
}

# Main script execution
# Main script execution
main() {
    PARSE_OUTPUT=$(parse_command_line_options "$@")
    if [ $? -eq 1 ]; then
        return 1
    fi
    read IS_DRY_RUN IS_SKIP_JULIA IS_SKIP_RSYNC SCRIPT_TO_RUN WORKSPACE_DIR DESTINATION <<< "$PARSE_OUTPUT"
    if [ $IS_SKIP_JULIA -eq 0 ]; then
        run_julia_script "$SCRIPT_TO_RUN"
        # ""
        awk '/nstep_out/ {print; print "    nstep_rdb        =  500"; next}1' output.dat > temp && mv temp output.dat
        awk '/min_cohes/ {print; print "    DII_ref          =  1e-18"; next}1' output.dat > temp && mv temp output.dat

    fi
    if [ $IS_SKIP_RSYNC -eq 0 ]; then
        ssh -t alaik@snellius.surf.nl "cd lam && cleanL.sh"
        perform_rsync_operation "$WORKSPACE_DIR" "$DESTINATION"
    fi
    ssh_and_run_script
}
main "$@"