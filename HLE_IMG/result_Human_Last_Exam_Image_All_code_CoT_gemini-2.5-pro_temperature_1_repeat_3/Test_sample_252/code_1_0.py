import cv2
import numpy as np
from PIL import Image
import requests
from io import BytesIO

# The URL of the original image provided in the problem
# This is the first image in the list of images attached.
image_url = "https://i.imgur.com/k2A9c97.png"

def load_image_from_url(url):
    """Loads an image from a URL into an OpenCV-compatible format."""
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        img = Image.open(BytesIO(response.content)).convert('RGB')
        return cv2.cvtColor(np.array(img), cv2.COLOR_RGB2BGR)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading image: {e}")
        return None

# Load the original image from the URL
original_image = load_image_from_url(image_url)

if original_image is not None:
    # Option C: Apply a Gaussian filter with a 7x7 pixel kernel.
    # The sigmaX and sigmaY parameters determine the standard deviation of the Gaussian kernel.
    # If set to 0, they are calculated from the kernel size.
    kernel_size = (7, 7)
    processed_image_c = cv2.GaussianBlur(original_image, kernel_size, 0)

    # Save the output image to a file so the user can view it.
    output_filename = 'processed_image_gaussian_7x7.png'
    cv2.imwrite(output_filename, processed_image_c)

    print(f"The original image was processed using a Gaussian filter with a {kernel_size[0]}x{kernel_size[1]} kernel.")
    print(f"The resulting image has been saved as '{output_filename}'.")
    print("Please compare this generated image with the 'Processed Image' in the problem description to verify the result.")
else:
    print("Could not load the image. Please check the URL or your internet connection.")
