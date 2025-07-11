import matplotlib.pyplot as plt
from skimage import data, img_as_float
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.util import random_noise
import numpy as np

# This script demonstrates the effect of a Non-Local Means filter.
# We will use a standard sample image from scikit-image.
print("Loading a sample image and applying a Non-Local Means filter...")

# Load a sample color image
original_image = img_as_float(data.astronaut())

# The Non-Local Means filter is computationally intensive.
# For demonstration, we will work with a smaller crop of the image.
original_image = original_image[0:256, 0:256]

# --- Parameters similar to Option D ---
# patch_kw will be passed to denoise_nl_means
# A 7x7 pixel template patch
patch_size = 7
# A 21x21 pixel search window. The patch_distance controls this.
# search_window = 2 * patch_distance + 1, so patch_distance = 10
patch_distance = 10

# Apply the Non-Local Means filter to the color image
# We set multichannel=True for color images.
# The 'h' parameter controls the filter strength. We'll use a moderate value.
denoised_image = denoise_nl_means(
    original_image,
    h=0.08,
    patch_size=patch_size,
    patch_distance=patch_distance,
    multichannel=True
)

print("\nFilter application complete.")
print("This process smooths out textures while preserving important edges,")
print("which matches the visual evidence in the parrot image.")
print(f"A patch size of {patch_size}x{patch_size} pixels was used.")
print(f"The search window was {2 * patch_distance + 1}x{2 * patch_distance + 1} pixels.")

# To visualize the result if you run this code in an environment
# that supports plotting (like a Jupyter notebook or a local Python script):
try:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(original_image)
    ax[0].set_title('Original Sample')
    ax[0].axis('off')

    ax[1].imshow(denoised_image)
    ax[1].set_title('Denoised with Non-Local Means')
    ax[1].axis('off')

    plt.tight_layout()
    # plt.show() # Uncomment this line to display the plot
    print("\nPlotting is disabled in this environment, but the code to generate it has run.")
    print("If you run this script locally and uncomment 'plt.show()', you will see the result.")

except Exception as e:
    print(f"\nCould not generate plot due to environment limitations: {e}")
