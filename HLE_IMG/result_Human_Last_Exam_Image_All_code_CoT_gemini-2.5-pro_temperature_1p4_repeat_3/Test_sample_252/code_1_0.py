import cv2
import numpy as np
from skimage import io
import matplotlib.pyplot as plt

# --- User Configuration ---
# You need to have an image named 'parrot.png' in the same directory.
# The code will generate 'processed_parrot.png' as the output.
IMAGE_URL = 'https://i.imgur.com/e44T8aC.png'
INPUT_FILENAME = 'parrot.png'
OUTPUT_FILENAME = 'processed_parrot.png'

# --- Main Script ---
try:
    # Download the image if it doesn't exist
    import os
    if not os.path.exists(INPUT_FILENAME):
        print(f"Downloading image from {IMAGE_URL} to {INPUT_FILENAME}...")
        original_image = io.imread(IMAGE_URL)
        # Convert RGBA to RGB if necessary
        if original_image.shape[2] == 4:
            original_image = original_image[:, :, :3]
        plt.imsave(INPUT_FILENAME, original_image)
        print("Download complete.")

    # Read the original image using OpenCV
    # OpenCV reads images in BGR format by default
    img = cv2.imread(INPUT_FILENAME)

    if img is None:
        print(f"Error: Could not read the image file '{INPUT_FILENAME}'.")
        print("Please make sure the image is in the correct directory.")
    else:
        print("Applying Non-Local Means filter...")
        # Parameters from option D:
        # h: Filter strength. A larger h removes more noise but can make the image look blurry.
        h = 10
        # templateWindowSize: Size of the patch used to compute weights. Must be odd. (7x7)
        templateWindowSize = 7
        # searchWindowSize: Size of the window where we search for similar patches. Must be odd. (21x21)
        searchWindowSize = 21

        # Apply the Non-Local Means filter for color images
        processed_img = cv2.fastNlMeansDenoisingColored(img, None, h, h, templateWindowSize, searchWindowSize)

        # Save the processed image
        cv2.imwrite(OUTPUT_FILENAME, processed_img)

        print(f"Processing complete. The processed image has been saved as '{OUTPUT_FILENAME}'.")

except ImportError:
    print("Error: Required libraries not found.")
    print("Please install them using: pip install opencv-python numpy scikit-image matplotlib")
except Exception as e:
    print(f"An error occurred: {e}")
