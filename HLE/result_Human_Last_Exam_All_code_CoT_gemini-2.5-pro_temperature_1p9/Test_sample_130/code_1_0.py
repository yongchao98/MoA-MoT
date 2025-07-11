import numpy as np
from scipy.ndimage import median_filter

def demonstrate_median_flow():
    """
    Simulates the repeated application of a local median filter on a binary image
    to show what happens to edges as t -> infinity.
    """

    # --- Parameters for the simulation ---
    # These define the specific instance of the process.
    # The 'equation' being modeled is: Image(t+1) = MedianFilter(Image(t), radius=delta)
    N = 20      # Image dimension (N x N)
    delta = 3   # Radius for the local median function (delta << N)
    max_iterations = 30 # Number of times to apply the operator (t)

    # --- Initial Image Setup ---
    # Create a binary image I in {0,1}^(N x N)
    # Here, we create a white square on a black background. The edges
    # of this square are the boundaries we will observe.
    image = np.zeros((N, N), dtype=np.uint8)
    square_start = (N - 10) // 2
    square_end = square_start + 10
    image[square_start:square_end, square_start:square_end] = 1

    print("--- Simulation of Repeated Local Median Filtering ---")
    print(f"This script demonstrates the behavior of the function: f_t+1(x) = Median_||x-y||<=d(f_t(y))")
    print(f"We will use the following numerical values for the parameters:")
    print(f"Image Size (N): {N}")
    print(f"Neighborhood Radius (delta): {delta}")
    print("-" * 50)
    
    # --- The Simulation Loop ---
    # We track the "area" of the white shape by summing all pixel values.
    # A decrease in area indicates the shape is shrinking.
    
    current_image = image.copy()
    initial_area = np.sum(current_image)
    print(f"Area of the white shape at t=0: {initial_area}")

    # The neighborhood for the median is a circle of radius delta.
    # We create a 'footprint' mask for scipy's median_filter to represent this circle.
    y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
    footprint = (x**2 + y**2 <= delta**2)

    for t in range(1, max_iterations + 1):
        # Apply the local median operator
        current_image = median_filter(current_image, footprint=footprint, mode='constant', cval=0)
        current_area = np.sum(current_image)
        
        print(f"Area of the white shape at t={t}: {current_area}")

        if current_area == 0:
            print("\nObservation: The area has reached zero.")
            print(f"The shape, and therefore its edges, has vanished completely after {t} iterations.")
            break
            
    if np.sum(current_image) > 0:
        print("\nObservation: The area is decreasing with each iteration.")
        print("The shape is shrinking. If we continued for more iterations, it would eventually vanish.")

    print("-" * 50)
    print("Final Conclusion for t -> infinity:")
    print("The repeated application of the local median filter causes edges to smooth out and convex shapes to shrink.")
    print("Ultimately, as t approaches infinity, any closed edges will disappear, and the image will become uniform.")

demonstrate_median_flow()