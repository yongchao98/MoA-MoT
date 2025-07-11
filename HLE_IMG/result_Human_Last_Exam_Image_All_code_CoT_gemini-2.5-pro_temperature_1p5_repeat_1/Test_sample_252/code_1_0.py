import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.util import img_as_ubyte, img_as_float

# --- Instructions for the user ---
# 1. Make sure you have the required libraries installed:
#    pip install scikit-image matplotlib imageio numpy
#
# 2. Download the original image from the problem description and save it
#    as 'parrot_original.png' in the same directory as this script.
#
# 3. If you cannot download the image, this code will generate a placeholder.

try:
    # Load the original image
    original_image = imageio.imread('parrot_original.png')
    # The provided image has an alpha channel, remove it for processing
    if original_image.shape[2] == 4:
        original_image = original_image[:, :, :3]
except FileNotFoundError:
    print("Image 'parrot_original.png' not found.")
    print("This script requires the original image file to run.")
    # As a fallback, create a noisy sample image to demonstrate the filter.
    # Note: The effect will be more clear on the actual parrot image.
    c = np.zeros([128, 128, 3])
    c[32:-32, 32:-32, :] = 1
    c = c + 0.3 * np.random.randn(*c.shape)
    c[c < 0] = 0
    c[c > 1] = 1
    original_image = img_as_ubyte(c)


# Convert image to floating point for processing
original_float = img_as_float(original_image)

# Parameters for the Non-Local Means filter, as specified in Option D
# patch_size: Size of patches used for comparison.
# patch_distance: Maximum distance to search for patches. Corresponds to search window size.
# Search window size = 2 * patch_distance + 1. So 21x21 window -> patch_distance = 10
patch_kw = dict(patch_size=7,      # 7x7 pixel template patch
                patch_distance=10, # Search window is 2*10+1 = 21x21
                multichannel=True)

# Estimate the noise standard deviation from the image.
# This helps in setting the filtering strength `h`.
# A higher `h` leads to more smoothing. We'll set it manually for a strong effect.
# sigma_est = np.mean(estimate_sigma(original_float, multichannel=True))
# The `h` parameter controls filter strength. It is often set relative to sigma.
# For a visible smoothing effect on a non-noisy image, we set a value manually.
h_value = 0.08 # Adjust this value to see different smoothing strengths

# Apply the Non-Local Means filter
print("Applying Non-Local Means filter...")
print(f"Template patch size: {patch_kw['patch_size']}x{patch_kw['patch_size']} pixels")
print(f"Search window size: {2 * patch_kw['patch_distance'] + 1}x{2 * patch_kw['patch_distance'] + 1} pixels")

denoised_image = denoise_nl_means(original_float, h=h_value, fast_mode=True, **patch_kw)

# Display the results
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(original_image)
ax[0].set_title('Original Image')
ax[0].axis('off')

ax[1].imshow(denoised_image)
ax[1].set_title('Processed with Non-Local Means')
ax[1].axis('off')

plt.tight_layout()
plt.show()
