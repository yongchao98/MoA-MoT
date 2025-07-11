import cv2
import numpy as np
import requests
from PIL import Image
import io

def process_image():
    """
    This function downloads an image of a parrot, applies a Gaussian blur
    with a 7x7 kernel, and saves the result. This demonstrates the effect
    described in Option C.
    """
    try:
        # URL of the original image
        image_url = "https://i.imgur.com/k91J0nQ.png"
        
        # Download the image
        print("Downloading the original image...")
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Read the image from the response content
        pil_image = Image.open(io.BytesIO(response.content)).convert('RGB')
        
        # Convert PIL image to OpenCV format (NumPy array)
        # PIL uses RGB, OpenCV uses BGR, so we need to convert color channels
        original_image = cv.cvtColor(np.array(pil_image), cv.COLOR_RGB2BGR)
        
        # Define the kernel size for the Gaussian filter
        kernel_size = (7, 7)
        
        # Apply the Gaussian filter
        # The third argument is sigmaX. When set to 0, it's calculated from the kernel size.
        print(f"Applying a Gaussian filter with a {kernel_size[0]}x{kernel_size[1]} kernel...")
        processed_image = cv.GaussianBlur(original_image, kernel_size, 0)
        
        # Save the processed image to a file
        output_filename = "processed_parrot.png"
        cv.imwrite(output_filename, processed_image)
        
        print(f"Processing complete. The processed image has been saved as '{output_filename}'.")

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    process_image()
