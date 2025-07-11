import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import data

def demonstrate_gaussian_blur():
    """
    This function demonstrates the effect of a Gaussian filter on an image,
    which is the most likely operation performed on the parrot image.
    """
    # Load a sample image with good texture and detail
    # The 'chelsea' image (a cat) from scikit-image is a standard example.
    # OpenCV loads images in BGR order, but matplotlib expects RGB.
    # We will convert it for correct color display.
    original_image_bgr = data.chelsea()
    original_image = cv2.cvtColor(original_image_bgr, cv2.COLOR_BGR2RGB)

    # Define the kernel size for the Gaussian filter.
    # The problem specifies a 7x7 kernel.
    kernel_width = 7
    kernel_height = 7
    kernel_size = (kernel_width, kernel_height)

    # Apply the Gaussian filter.
    # The third argument (sigmaX) is the standard deviation.
    # Setting it to 0 lets OpenCV calculate it automatically from the kernel size.
    print(f"Applying Gaussian filter with a {kernel_width}x{kernel_height} kernel...")
    processed_image = cv2.GaussianBlur(original_image, kernel_size, 0)

    # Display the original and processed images side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].imshow(original_image)
    axes[0].set_title('Original Image')
    axes[0].axis('off')

    axes[1].imshow(processed_image)
    axes[1].set_title(f'Processed with {kernel_width}x{kernel_height} Gaussian Filter')
    axes[1].axis('off')

    plt.tight_layout()
    plt.show()

demonstrate_gaussian_blur()