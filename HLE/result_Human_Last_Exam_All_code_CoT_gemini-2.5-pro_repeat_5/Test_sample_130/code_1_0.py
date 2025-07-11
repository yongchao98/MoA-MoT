import numpy as np
from scipy.ndimage import median_filter
from skimage.morphology import disk

# --- Parameters ---
# Image size
N = 15
# Radius of the median filter neighborhood
delta = 2

# --- 1. Create the initial image (t=0) ---
# A white square on a black background
image = np.zeros((N, N), dtype=int)
image[4:11, 4:11] = 1

print("--- Image at t=0 ---")
print("The initial image is a square with sharp edges. Each number represents a pixel:")
for row in image:
    # The 'end' and 'flush' parameters are used for compatibility across different terminals.
    print(' '.join(map(str, row)), end='\n', flush=True)

# --- 2. Define the filter ---
# Create a circular neighborhood (footprint) for the filter, as specified by ||x-y|| <= delta
footprint = disk(delta)

# --- 3. Apply the filter iteratively ---
# Apply the filter for the first time
image_t1 = median_filter(image, footprint=footprint, mode='constant', cval=0)

print("\n--- Image at t=1 ---")
print("After one iteration, the sharp corners of the square begin to round off:")
for row in image_t1:
    print(' '.join(map(str, row)), end='\n', flush=True)

# Apply the filter a few more times to see more significant changes
image_t5 = image.copy()
for _ in range(5):
    image_t5 = median_filter(image_t5, footprint=footprint, mode='constant', cval=0)

print("\n--- Image at t=5 ---")
print("After five iterations, the shape is visibly smoother and more circular:")
for row in image_t5:
    print(' '.join(map(str, row)), end='\n', flush=True)
