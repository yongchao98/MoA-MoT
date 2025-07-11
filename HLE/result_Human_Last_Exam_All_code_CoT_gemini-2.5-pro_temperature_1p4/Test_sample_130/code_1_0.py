import numpy as np
from scipy.signal import medfilt2d

def calculate_perimeter(image):
    """
    Calculates the total length of edges in a binary image by summing the
    differences between adjacent pixels.
    """
    # Calculate differences between horizontally adjacent pixels
    h_diff = np.abs(image[:, 1:] - image[:, :-1])
    # Calculate differences between vertically adjacent pixels
    v_diff = np.abs(image[1:, :] - image[:-1, :])
    # Total perimeter is the sum of all these differences
    return np.sum(h_diff) + np.sum(v_diff)

# 1. Define image parameters for the simulation
N = 100  # Image is N x N pixels
delta = 2  # The radius for the median filter
kernel_size = 2 * delta + 1  # The filter looks at a (2*delta+1) x (2*delta+1) square
num_iterations = 30 # Number of times to apply the filter

# 2. Create an initial image: a 50x50 white square on a black background
image = np.zeros((N, N), dtype=np.uint8)
start, end = 25, 75
image[start:end, start:end] = 1

print("--- Simulation of repeated median filtering on a binary image ---")
print(f"Initial shape: A {end-start}x{end-start} square on a {N}x{N} grid.")
print(f"Filter kernel size: {kernel_size}x{kernel_size}")
print("-" * 60)
print("t     | Perimeter | Area (White Pixels)")
print("-" * 60)

# 3. Run the simulation loop
current_image = image.copy()
for t in range(num_iterations):
    perimeter = calculate_perimeter(current_image)
    white_pixels = np.sum(current_image)

    # The area of a square is (end-start)^2 = 50^2 = 2500
    # The perimeter of a square is 4 * (end-start) = 4 * 50 = 200
    print(f"{t:<5d} | {perimeter:<9.0f} | {white_pixels:<18d}")

    # If the shape has disappeared, stop the simulation
    if white_pixels == 0:
        print("\nShape has disappeared. Halting simulation.")
        break

    # Apply the median filter to get the image for the next step
    # The scipy function requires a float or int array, not boolean
    current_image = medfilt2d(current_image.astype(float), kernel_size=kernel_size).astype(np.uint8)
