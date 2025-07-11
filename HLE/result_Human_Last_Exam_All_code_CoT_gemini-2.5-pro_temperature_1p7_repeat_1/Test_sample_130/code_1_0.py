import numpy as np
from scipy.ndimage import median_filter

def simulate_median_filtering():
    """
    Simulates the repeated application of a local median filter on a binary image
    to observe the effect on edges over time.
    """
    # 1. Create a sample binary image (a 12x12 white square on a 30x40 black background)
    image_height, image_width = 30, 40
    img = np.zeros((image_height, image_width), dtype=np.int8)
    
    # Place a white square in the image
    square_size = 12
    start_row = (image_height - square_size) // 2
    start_col = (image_width - square_size) // 2
    img[start_row:start_row + square_size, start_col:start_col + square_size] = 1

    # 2. Define the neighborhood for the median filter.
    # The neighborhood is a disk of radius delta, as specified by ||x-y|| <= delta.
    delta = 3
    y, x = np.ogrid[-delta:delta + 1, -delta:delta + 1]
    # Create a circular footprint (structuring element)
    disk_footprint = (x**2 + y**2 <= delta**2)

    def print_image_and_stats(image, time_step):
        """Helper function to print the image and the number of white pixels."""
        print(f"\n--- Image at t={time_step} ---")
        # Print a character representation of the image
        for row in image:
            print(" ".join(['■' if p == 1 else '·' for p in row]))
        
        # Calculate and print the number of white pixels to quantify the change
        num_white_pixels = np.sum(image)
        print(f"Number of white pixels: {num_white_pixels}")
        if time_step == 0:
            print(f"(Initial area is {square_size}x{square_size} = {square_size**2})")

    # --- Simulation ---
    
    # Print initial state
    print_image_and_stats(img, 0)
    
    # Apply the filter iteratively and show the result at different time steps
    time_steps_to_show = [5, 10, 15, 20]
    
    current_image = img.copy()
    for t in range(1, max(time_steps_to_show) + 1):
        # Apply the median filter. On a binary image, this acts as a majority filter.
        # mode='constant' with cval=0 handles pixels near the border by padding with black pixels.
        current_image = median_filter(current_image, footprint=disk_footprint, mode='constant', cval=0)
        
        if t in time_steps_to_show:
            print_image_and_stats(current_image, t)
            # Stop if the object has completely disappeared
            if np.sum(current_image) == 0:
                print("\nThe white square has completely disappeared.")
                break

if __name__ == '__main__':
    simulate_median_filtering()