import cv2
import numpy as np

# --- User Action Required ---
# Please save the original image from the problem description as 'original_parrot.png'
# and place it in the same folder as this script before running.

def apply_non_local_means_filter():
    """
    Loads an image, applies a Non-Local Means filter, and saves the result.
    """
    input_filename = 'original_parrot.png'
    output_filename = 'processed_parrot.png'

    # Step 1: Load the original image
    try:
        original_image = cv2.imread(input_filename)
        if original_image is None:
            print(f"Error: Could not load the image '{input_filename}'.")
            print("Please make sure the file exists and is in the correct directory.")
            return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    print("Successfully loaded original image.")

    # Step 2: Define filter parameters based on Option D
    # Parameters for the cv2.fastNlMeansDenoisingColored function:
    # h and hColor: Filter strength. 10 is a good value for visible smoothing.
    denoising_strength = 10
    # templateWindowSize: Size of the patch used to compute weights. Must be odd.
    template_patch_size = 7
    # searchWindowSize: Size of the window where patches are searched. Must be odd.
    search_window_size = 21

    print("\nApplying Non-Local Means filter with the following parameters:")
    print(f"Template Patch Size  = {template_patch_size}x{template_patch_size} pixels")
    print(f"Search Window Size   = {search_window_size}x{search_window_size} pixels")
    print(f"Denoising Strength   = {denoising_strength}")

    # Step 3: Apply the filter
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image,
        None,
        denoising_strength,
        denoising_strength,
        template_patch_size,
        search_window_size
    )

    # Step 4: Save the processed image
    try:
        cv2.imwrite(output_filename, processed_image)
        print(f"\nProcessing complete. The resulting image has been saved as '{output_filename}'.")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

if __name__ == '__main__':
    # To run this code, you need to have opencv-python installed:
    # pip install opencv-python
    apply_non_local_means_filter()