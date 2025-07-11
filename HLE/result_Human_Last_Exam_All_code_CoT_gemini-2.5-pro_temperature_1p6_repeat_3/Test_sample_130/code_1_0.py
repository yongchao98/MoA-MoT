import numpy as np
from scipy.ndimage import median_filter

def create_plus_sign_image(size=20, arm_width=6):
    """Creates a binary image with a plus sign shape."""
    image = np.zeros((size, size), dtype=int)
    center = size // 2
    half_width = arm_width // 2
    # Define boundaries for the arms to make a plus sign
    arm_start = size // 5
    arm_end = 4 * size // 5
    # Horizontal arm
    image[center - half_width : center + half_width, arm_start:arm_end] = 1
    # Vertical arm
    image[arm_start:arm_end, center - half_width : center + half_width] = 1
    return image

def print_image_ascii(image, title):
    """Prints the image array as ASCII art for clear visualization."""
    print(f"--- {title} ---")
    # Using '■' for 1 and a space for 0 to represent the shape.
    for row in image:
        print("".join(['■ ' if pixel == 1 else '  ' for pixel in row]))
    print("\n")

# --- Parameters ---
image_size = 20
delta = 2           # Radius of the median filter neighborhood
total_iterations = 10 # Number of times the filter is applied

# --- Setup ---
# Create the initial image with a plus sign
image = create_plus_sign_image(size=image_size)

# Create the disk-shaped structuring element for the filter's neighborhood
y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
disk_kernel = x**2 + y**2 <= delta**2

# --- Simulation & Output ---
print("This simulation demonstrates the effect of an iterative median filter on a binary image.")
print("The initial shape is a plus sign ('■').\n")

# Print the initial state of the image
print_image_ascii(image, "Image at t = 0")

# Apply the median filter iteratively
processed_image = np.copy(image)
for t in range(total_iterations):
    processed_image = median_filter(processed_image, footprint=disk_kernel)

# Print the final state of the image after all iterations
print_image_ascii(processed_image, f"Image at t = {total_iterations} (Final State Visualized)")

print(f"--- Final Image Numerical Array (t = {total_iterations}) ---")
# As requested to "output each number in the final equation!",
# the final numerical state of the image is printed below.
for row in processed_image:
    print(' '.join(map(str, row)))