import cv2
import numpy as np
import requests
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt

def demonstrate_gaussian_blur():
    """
    Downloads the parrot image, applies a Gaussian blur as described in Option C,
    and displays the original and processed images for comparison.
    """
    try:
        # The URL points to the first cropped image of the parrot.
        image_url = 'https://i.imgur.com/zV8R9X1.png'
        print(f"Downloading image from: {image_url}")
        response = requests.get(image_url)
        # Raise an error if the download fails
        response.raise_for_status()

        # Read the image data into a Pillow Image and then convert to an OpenCV format (NumPy array)
        img_pil = Image.open(BytesIO(response.content)).convert('RGB')
        original_image = np.array(img_pil)
        # OpenCV processes images in BGR color order, so we need to convert from RGB
        original_image_bgr = cv2.cvtColor(original_image, cv2.COLOR_RGB2BGR)

        # Define the kernel size for the Gaussian filter
        kernel_width = 7
        kernel_height = 7
        kernel_size = (kernel_width, kernel_height)
        
        print(f"Applying Gaussian filter with a {kernel_width}x{kernel_height} kernel...")

        # Apply the Gaussian filter. A sigma of 0 tells OpenCV to calculate it from the kernel size.
        processed_image_bgr = cv2.GaussianBlur(original_image_bgr, kernel_size, 0)

        # Convert the processed image back to RGB for display with Matplotlib
        processed_image_rgb = cv2.cvtColor(processed_image_bgr, cv2.COLOR_BGR2RGB)

        # Display the results
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        
        axes[0].imshow(original_image)
        axes[0].set_title("Original Image")
        axes[0].axis('off')

        axes[1].imshow(processed_image_rgb)
        axes[1].set_title(f"Processed: {kernel_width}x{kernel_height} Gaussian Blur")
        axes[1].axis('off')

        plt.tight_layout()
        plt.show()

    except requests.exceptions.RequestException as e:
        print(f"Failed to download image. Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the demonstration
demonstrate_gaussian_blur()