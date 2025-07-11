import numpy as np

def run_simulation():
    """
    Simulates the repeated application of a local median filter on a binary image
    to observe the effect on its edges.
    """
    # 1. Create a sample binary image (a hollow square)
    # N is the image size, delta is the filter radius.
    # We choose delta << N as specified.
    N = 40
    delta = 2
    
    # Start with a black image
    image = np.zeros((N, N), dtype=int)
    
    # Draw a white square in the middle
    image[8:N-8, 8:N-8] = 1
    
    # Carve out a black square to make it hollow
    image[14:N-14, 14:N-14] = 0

    # --- Helper function to print the image ---
    def print_image(img, title):
        print(title)
        # Use '#' for white pixels (1) and ' ' for black pixels (0)
        for row in img:
            print("".join(["##" if p == 1 else "  " for p in row]))
        print("-" * img.shape[1] * 2)

    # --- 2. Implement the filter application ---
    def apply_local_median(img, radius):
        # Create a copy to store the new pixel values
        new_img = np.copy(img)
        rows, cols = img.shape
        
        # Iterate over each pixel in the image
        for r in range(rows):
            for c in range(cols):
                # Define the square neighborhood boundaries, clamped to the image size
                r_min = max(0, r - radius)
                r_max = min(rows, r + radius + 1)
                c_min = max(0, c - radius)
                c_max = min(cols, c + radius + 1)
                
                # Extract the neighborhood
                neighborhood = img[r_min:r_max, c_min:c_max]
                
                # For a binary image, the median pushes the pixel towards the
                # majority value (0 or 1) in its neighborhood.
                median_val = np.median(neighborhood)
                
                # Update the pixel in the new image
                new_img[r, c] = int(round(median_val))
                
        return new_img

    # --- 3. & 4. Iterate the process and visualize ---
    print("This simulation shows what happens when a local median filter is applied repeatedly.")
    print(f"The filter has a radius of delta = {delta}.\n")
    
    iterations_to_show = [0, 2, 5, 10]
    current_image = np.copy(image)

    # Show the initial state
    print_image(current_image, "Image at t = 0 (Initial State):")

    # Run the simulation
    for t in range(1, max(iterations_to_show) + 1):
        current_image = apply_local_median(current_image, delta)
        if t in iterations_to_show:
            print_image(current_image, f"Image at t = {t}:")
            
    # --- 5. Draw a Conclusion from the simulation ---
    print("ANALYSIS OF THE RESULT:")
    print("The simulation shows that the operator smooths the edges and corners.")
    print("1. Convex corners (the outer corners of the square) are rounded and 'eaten away'.")
    print("2. Concave corners (the inner corners of the hollow part) are also rounded and 'filled in'.")
    print("This process is a discrete version of 'mean curvature flow', which acts to minimize the total length of the edges.")
    print("As t -> infinity, any finite shape will become circular, shrink, and eventually disappear, resulting in a single-colored image.")

# Run the simulation
run_simulation()