import os
import pandas as pd
import urllib.request
from tqdm import tqdm

# --- CONFIGURATION ---
# 1. The name of your Runinfo file.
RUNINFO_FILE = "RunInfo.csv" # The file you downloaded

# 2. The corrected column names from your file.
COLUMN_RUN = "Run"
COLUMN_URL = "Download_path"      # CORRECTED
COLUMN_SAMPLE = "SampleName"        # CORRECTED

# 3. The name of the output directory.
OUTPUT_DIR = "raw_fastq"
# --- END OF CONFIGURATION ---


class TqdmUpTo(tqdm):
    """Provides a progress bar hook for urllib.request.urlretrieve."""
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_data():
    """
    Reads the Runinfo file and downloads paired-end FASTQ files.
    """
    if not os.path.exists(RUNINFO_FILE):
        print(f"Error: Runinfo file '{RUNINFO_FILE}' not found.")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Reading metadata from: {RUNINFO_FILE}")
    # CORRECTED: Changed separator to a comma ','
    df = pd.read_csv(RUNINFO_FILE, sep=',')

    print(f"Found {len(df)} samples to process. Starting downloads...")
    print("-----------------------------------------------------")

    for index, row in df.iterrows():
        sample_name = row[COLUMN_SAMPLE]
        run_accession = row[COLUMN_RUN]
        fastq_urls = row[COLUMN_URL]

        print(f"Processing Sample: {sample_name} (Run: {run_accession})")

        try:
            # CORRECTED: Split the URLs string by the pipe '|' character
            url_r1, url_r2 = fastq_urls.split('|')
        except (ValueError, AttributeError):
            print(f"  -> Warning: Could not parse URLs for {sample_name}. Skipping.")
            print("-----------------------------------------------------")
            continue

        filepath_r1 = os.path.join(OUTPUT_DIR, f"{sample_name}_R1.fastq.gz")
        filepath_r2 = os.path.join(OUTPUT_DIR, f"{sample_name}_R2.fastq.gz")

        for url, filepath in [(url_r1, filepath_r1), (url_r2, filepath_r2)]:
            filename = os.path.basename(filepath)
            if not os.path.exists(filepath):
                print(f"  -> Downloading {filename}...")
                try:
                    with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=filename) as t:
                        urllib.request.urlretrieve(url, filename=filepath, reporthook=t.update_to)
                except Exception as e:
                    print(f"\n  -> Error downloading {url}: {e}")
                    if os.path.exists(filepath):
                        os.remove(filepath)
            else:
                print(f"  -> File {filename} already exists. Skipping.")
        
        print("-----------------------------------------------------")

    print("✅ All downloads are complete.")


if __name__ == "__main__":
    download_data()