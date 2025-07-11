import numpy as np

def apply_local_median(image, delta):
    """
    Applies a local median filter to a binary image.
    
    For each pixel in the input image, it considers a square neighborhood of
    size (2*delta + 1) x (2*delta + 1) centered at that pixel. The new
    value of the pixel is determined by the median of the values in its
    neighborhood. For a binary {0, 1} image, this is equivalent to a
    majority vote: if more than half the pixels in the neighborhood are 1,
    the output pixel is 1; otherwise, it is 0.

    Args:
        image (np.ndarray): The input binary image.
        delta (int): The radius of the square neighborhood. The window
                     size will be (2*delta + 1).

    Returns:
        np.ndarray: The image after applying the median filter.
    """
    h, w = image.shape
    # Pad the image with zeros to handle the borders correctly.
    padded_image = np.pad(image, delta, mode='constant', constant_values=0)
    new_image = np.zeros_like(image)

    # Iterate over each pixel of the original image to compute its new value.
    for i in range(h):
        for j in range(w):
            # Define the slice for the neighborhood in the padded image.
            # The neighborhood for pixel (i, j) is at a corresponding
            # location in the padded image.
            neighborhood = padded_image[i : i + 2 * delta + 1,
                                        j : j + 2 * delta + 1]
            
            # The median in a binary array is 1 if the sum is greater
            # than half the size of the array, and 0 otherwise.
            if np.sum(neighborhood) > (neighborhood.size / 2.0):
                new_image[i, j] = 1
            else:
                new_image[i, j] = 0
                
    return new_image

# --- Simulation Parameters ---
N = 20  # Image is N x N
delta = 1 # Neighborhood radius. delta << N
iterations_to_show = [0, 1, 2, 3, 5, 8]

# --- Initial Image Setup ---
# Create an initial binary image: a white square on a black background.
image = np.zeros((N, N), dtype=int)
square_size = 10
start_idx = (N - square_size) // 2
end_idx = start_idx + square_size
image[start_idx:end_idx, start_idx:end_idx] = 1

# --- Simulation and Explanation ---
print("This script simulates the iterative application of a local median filter on a binary image.")
print("We start with a white square to observe the effect on its edges and corners.\n")

current_image = image.copy()
max_iter = max(iterations_to_show) + 1 if iterations_to_show else 1

for t in range(max_iter):
    if t in iterations_to_show:
        print("-" * (2 * N + 5))
        print(f"Image at iteration t={t}:")
        
        # Check if the image has become all black
        if np.sum(current_image) == 0 and t > 0:
            print("(The image has become completely black)")
        
        # Print the image using '##' for white and '..' for black for better visibility
        for row in current_image:
            print(" ".join(['##' if p == 1 else '..' for p in row]))
        
        if np.sum(current_image) == 0 and t > 0:
             break
    
    current_image = apply_local_median(current_image, delta)


print("-" * (2 * N + 5))
print("\n--- Final Analysis ---")
print("The simulation shows how the edges of the white square evolve over time:")
print("1. Smoothing of Corners: The sharp 90-degree corners of the square are rounded off first. This is because pixels at sharp convex corners have more neighbors outside the shape than inside, so the majority vote turns them black.")
print("2. Contraction of Edges: After the corners are smoothed, the straight edges begin to move inwards, causing the entire shape to shrink.")
print("3. Annihilation: As t -> infinity, this process continues. The shape gets progressively smaller until it eventually collapses into a point and disappears entirely. The final image becomes uniform (all black).")
print("\nThis behavior is a discrete approximation of a geometric process called 'motion by mean curvature'. The edges move with a velocity proportional to their curvature, causing them to smooth out, contract, and ultimately vanish.")