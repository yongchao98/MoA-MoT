import numpy as np

def run_median_filter_simulation():
    """
    Simulates the iterative application of a local median filter on a binary image
    and prints the evolution of the shape's properties.
    """
    # 1. SETUP
    N = 100  # Image size (N x N)
    delta = 5  # Radius of the median filter neighborhood
    
    # Create an initial image with a white square on a black background
    image = np.zeros((N, N), dtype=np.uint8)
    start, end = 25, 75
    image[start:end, start:end] = 1

    # Pre-calculate a grid of coordinates for distance calculation
    coords = np.arange(N)
    xx, yy = np.meshgrid(coords, coords)

    def apply_local_median(img, radius):
        """Applies one step of the local median filter."""
        new_img = np.copy(img)
        radius_sq = radius**2
        for r in range(N):
            for c in range(N):
                # Define the circular neighborhood
                dist_sq = (xx - c)**2 + (yy - r)**2
                mask = dist_sq <= radius_sq
                
                # Get neighbor values
                neighbors = img[mask]
                
                # Calculate the median for binary data (majority rule)
                # np.median returns 0.5 for a tie; we'll map that to 0.
                median_val = np.median(neighbors)
                new_img[r, c] = 1 if median_val > 0.5 else 0
        return new_img

    def calculate_properties(img):
        """Calculates the area and perimeter of the white shape."""
        area = np.sum(img)
        
        # Calculate perimeter by finding white pixels adjacent to black pixels
        perimeter = 0
        # Padded version to handle boundaries easily
        padded_img = np.pad(img, pad_width=1, mode='constant', constant_values=0)
        for r_pad in range(1, N + 1):
            for c_pad in range(1, N + 1):
                if padded_img[r_pad, c_pad] == 1: # If the pixel is white
                    # Check 4-connectivity neighbors in the padded image
                    if (padded_img[r_pad-1, c_pad] == 0 or
                        padded_img[r_pad+1, c_pad] == 0 or
                        padded_img[r_pad, c_pad-1] == 0 or
                        padded_img[r_pad, c_pad+1] == 0):
                        perimeter += 1
        return area, perimeter

    # 2. SIMULATION
    print("--- Starting Simulation ---")
    
    current_image = image
    timesteps_to_print = [0, 1, 2, 5, 10, 20, 30]

    for t in range(max(timesteps_to_print) + 1):
        if t in timesteps_to_print:
            area, perimeter = calculate_properties(current_image)
            if area == 0:
                print(f"t={t}: Shape has vanished. Area=0, Perimeter=0")
                break
            print(f"t={t}: Area={area}, Perimeter={perimeter}")

        # Apply the filter for the next step
        next_image = apply_local_median(current_image, delta)

        # If the image stops changing, break the loop
        if np.array_equal(current_image, next_image):
            print(f"\nImage has stabilized at t={t+1}.")
            break
        
        current_image = next_image
        
    print("--- Simulation Finished ---")

# Run the simulation
run_median_filter_simulation()
