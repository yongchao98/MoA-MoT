import numpy as np
from scipy.ndimage import median_filter

def print_image(img, title):
    """Prints a text representation of the binary image and its white pixel count."""
    print(title)
    for row in img:
        print(" ".join(['##' if pixel == 1 else '..' for pixel in row]))
    print(f"White Pixels: {np.sum(img)}\n")

# --- Simulation Parameters ---
IMAGE_SIZE = 24
RADIUS = 2 # The delta parameter (δ) for the median filter neighborhood

# --- 1. Create the Initial Image ---
# Create a black canvas
image = np.zeros((IMAGE_SIZE, IMAGE_SIZE), dtype=np.uint8)

# Draw a white "plus" sign in the middle
center = IMAGE_SIZE // 2
arm_width = 2
arm_length = 7
image[center - arm_length: center + arm_length, center - arm_width: center + arm_width] = 1
image[center - arm_width: center + arm_width, center - arm_length: center + arm_length] = 1

# --- 2. Create the Median Filter Footprint ---
# Create a circular neighborhood (footprint) for the filter
y, x = np.ogrid[-RADIUS:RADIUS + 1, -RADIUS:RADIUS + 1]
footprint = (x**2 + y**2 <= RADIUS**2).astype(np.uint8)

# --- 3. Run and Visualize the Simulation ---
# Define the time steps at which to print the image
timesteps_to_print = [0, 1, 3, 7, 12, 18]

current_image = image.copy()

print("--- Start of Median Filter Simulation ---")

for t in range(max(timesteps_to_print) + 1):
    if t in timesteps_to_print:
        print_image(current_image, f"Image at Iteration t = {t}")
    
    # Apply the local median filter
    current_image = median_filter(current_image, footprint=footprint, mode='constant', cval=0)

print("--- End of Simulation ---")
print("As you can see, the shape's sharp convex corners are rounded,")
print("its sharp concave corners are filled, and the entire shape")
print("shrinks over time, eventually disappearing.")
