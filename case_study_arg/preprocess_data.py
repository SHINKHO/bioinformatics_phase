import os
import sys
import glob
import subprocess

# --- CONFIGURATION ---
# Directory names
RAW_DIR = "raw_fastq"
CLEANED_DIR = "cleaned_fastq"

# Performance settings
THREADS = 8  # How many CPU threads to use

def find_conda_paths():
    """
    Automatically finds Trimmomatic paths within the active Conda environment.
    """
    # Get the active Conda environment's path
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        print("❌ Error: CONDA_PREFIX variable not found.", file=sys.stderr)
        print("      Please activate the Conda environment where Trimmomatic is installed.", file=sys.stderr)
        return None, None

    # Construct the expected paths based on standard Conda installation structure
    trimmomatic_jar = os.path.join(conda_prefix, "share", "trimmomatic", "trimmomatic.jar")
    adapters_fa = os.path.join(conda_prefix, "share", "trimmomatic", "adapters", "TruSeq3-PE-2.fa")

    # Check if the files actually exist at the constructed paths
    if not os.path.exists(trimmomatic_jar):
        print(f"❌ Error: Could not find 'trimmomatic.jar' in the Conda env.", file=sys.stderr)
        print(f"      Searched at: {trimmomatic_jar}", file=sys.stderr)
        print("      Please ensure Trimmomatic is installed (`conda install trimmomatic`).", file=sys.stderr)
        return None, None
    
    return trimmomatic_jar, adapters_fa


def run_command(command):
    """Runs a command using subprocess and checks for errors."""
    print(f"🚀 Executing: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error executing command. Aborting.", file=sys.stderr)
        print(f"STDERR:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Finds FASTQ files and runs Trimmomatic to clean them.
    """
    # --- Step 1: Automatically find Conda paths ---
    trimmomatic_jar_path, adapters_path = find_conda_paths()
    if not trimmomatic_jar_path:
        sys.exit(1) # Exit if paths couldn't be found
        
    print(f"✅ Found Trimmomatic JAR and adapter files in your Conda environment.")
    
    # --- Step 2: Set up directories and find files ---
    os.makedirs(CLEANED_DIR, exist_ok=True)
    r1_files = sorted(glob.glob(os.path.join(RAW_DIR, "*_R1.fastq.gz")))
    
    if not r1_files:
        print(f"❌ Error: No raw FASTQ files found in '{RAW_DIR}'.")
        return

    print(f"Found {len(r1_files)} samples to pre-process.")
    print("-----------------------------------------------------")

    # --- Step 3: Loop through samples and run Trimmomatic ---
    for r1_path in r1_files:
        base_name = os.path.basename(r1_path).replace("_R1.fastq.gz", "")
        print(f"Processing sample: {base_name}")
        
        r2_path = r1_path.replace("_R1.fastq.gz", "_R2.fastq.gz")
        
        r1_paired_out = os.path.join(CLEANED_DIR, f"{base_name}_R1_paired.fastq.gz")
        r1_unpaired_out = os.path.join(CLEANED_DIR, f"{base_name}_R1_unpaired.fastq.gz")
        r2_paired_out = os.path.join(CLEANED_DIR, f"{base_name}_R2_paired.fastq.gz")
        r2_unpaired_out = os.path.join(CLEANED_DIR, f"{base_name}_R2_unpaired.fastq.gz")

        trimmomatic_command = [
            "java", "-jar", trimmomatic_jar_path,
            "PE", "-threads", str(THREADS),
            r1_path, r2_path,
            r1_paired_out, r1_unpaired_out,
            r2_paired_out, r2_unpaired_out,
            f"ILLUMINACLIP:{adapters_path}:2:30:10",
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
        ]
        
        run_command(trimmomatic_command)
        print(f"✅ Finished cleaning {base_name}")
        print("-----------------------------------------------------")

    print(f"🎉 All samples have been processed. Cleaned files are in '{CLEANED_DIR}'.")


if __name__ == "__main__":
    main()