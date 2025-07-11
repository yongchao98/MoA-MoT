import numpy as np
import nibabel as nib
from io import BytesIO
import requests
import warnings

# Suppress a harmless warning from nibabel about file objects
warnings.filterwarnings("ignore", category=UserWarning, module='nibabel.filebasedimages')

def solve_neuro_query():
    """
    This script identifies the dmPFC sub-region most purely activated by emotion
    by querying the NeuroSynth meta-analytic database.
    It downloads the reverse inference Z-score map for the term 'emotion' and
    extracts the score at the coordinates of four dmPFC parcellations.
    The region with the highest Z-score is considered the most selectively activated.

    Required libraries:
    - numpy
    - nibabel
    - requests
    You can install them using: pip install numpy nibabel requests
    """
    try:
        # 1. Download the reverse inference map for 'emotion' from NeuroSynth
        print("Downloading NeuroSynth reverse inference map for 'emotion'...")
        url = 'https://github.com/neurosynth/neurosynth-data/raw/master/current_release/reverse_inference/emotion_z_rev.nii.gz'
        response = requests.get(url, timeout=20)
        response.raise_for_status()  # Ensure the download was successful
        print("Download complete.")

        # Load the Nifti image from the downloaded content in memory
        with BytesIO(response.content) as f:
            f.seek(0)
            img = nib.load(fileobj=f)

        # 2. Define representative MNI coordinates for the ROIs
        # These are based on well-known functional-anatomical distinctions in the literature
        rois = {
            'caudal-right': np.array([8, -12, 48]),
            'rostroventral': np.array([4, 34, 10]),
            'rostrodorsal': np.array([6, 26, 40]),
            'caudal-left': np.array([-6, -10, 50]),
        }

        # 3. Get the image data and the affine matrix to convert MNI to voxel space
        img_data = img.get_fdata()
        inv_affine = np.linalg.inv(img.affine)

        results = {}

        print("\nExtracting Z-scores at each region's coordinates...")
        # 4. Loop through ROIs, convert MNI to voxel coords, and get Z-score
        for name, mni_coords in rois.items():
            # Add a 1 to the end of the coords vector for affine transformation
            mni_coords_h = np.append(mni_coords, 1)
            
            # Apply the inverse affine to get voxel coordinates
            voxel_coords = np.round(inv_affine @ mni_coords_h).astype(int)
            x, y, z, _ = voxel_coords

            # Extract the Z-score from the image data at the voxel coordinate
            if 0 <= x < img_data.shape[0] and 0 <= y < img_data.shape[1] and 0 <= z < img_data.shape[2]:
                z_score = img_data[x, y, z]
                results[name] = z_score
            else:
                results[name] = np.nan # Mark as invalid if coordinate is outside the map

        # 5. Print the results and find the region with the highest score
        print("\n--- Results ---")
        print("Reverse inference Z-scores for the term 'emotion':")
        for name, score in results.items():
            print(f"Region: {name:<15} MNI: ({', '.join(map(str, rois[name]))}) ---> Z-score: {score:.4f}")

        # Find the region with the highest Z-score
        if not results:
             print("\nCould not retrieve any results.")
        else:
            # Filter out any potential NaN values before finding the max
            valid_results = {k: v for k, v in results.items() if not np.isnan(v)}
            if not valid_results:
                print("\nAll coordinates were outside the bounds of the map.")
            else:
                most_pure_region = max(valid_results, key=valid_results.get)
                max_score = valid_results[most_pure_region]
                print("\nConclusion:")
                print(f"The sub-region most purely (selectively) activated by emotion processing is '{most_pure_region}',")
                print(f"as it has the highest reverse inference Z-score of {max_score:.4f}.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred during download: {e}")
        print("Please check your internet connection and try again.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    solve_neuro_query()