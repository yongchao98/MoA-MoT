import math
import random

def calculate_area_of_R():
    """
    Estimates the area of the region R in the complex plane defined by
    floor(|z|) = |floor(z)| within the square from 0+0i to 6+6i.
    This is done using a Monte Carlo simulation.
    """
    
    # Number of random points to sample. A larger number gives a more accurate result.
    num_samples = 20_000_000
    
    # Counter for points falling inside the region R
    points_in_region = 0
    
    # The boundaries of the square in the complex plane
    x_min, x_max = 0.0, 6.0
    y_min, y_max = 0.0, 6.0
    
    # The simulation loop
    for _ in range(num_samples):
        # Generate a random complex number z = x + yi within the square
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)
        
        # Calculate floor(z) = floor(x) + i * floor(y)
        nx = math.floor(x)
        ny = math.floor(y)
        
        # Calculate |floor(z)|^2 = floor(x)^2 + floor(y)^2
        mod_floor_z_sq = nx * nx + ny * ny
        
        # Take the square root to get |floor(z)|
        mod_floor_z = math.sqrt(mod_floor_z_sq)
        
        # For the main equation to hold, |floor(z)| must be an integer.
        # We check this using a small tolerance for floating-point inaccuracies.
        if abs(mod_floor_z - round(mod_floor_z)) < 1e-9:
            # If it is an integer, let k be that integer value.
            k = round(mod_floor_z)
            
            # Calculate |z| = sqrt(x^2 + y^2)
            mod_z = math.sqrt(x * x + y * y)
            
            # Calculate floor(|z|)
            floor_mod_z = math.floor(mod_z)
            
            # Check the condition that defines the region R
            if floor_mod_z == k:
                points_in_region += 1
                
    # The area of the bounding square
    total_area = (x_max - x_min) * (y_max - y_min)
    
    # Estimate the area of R using the Monte Carlo formula
    estimated_area = (points_in_region / num_samples) * total_area
    
    # Output the results of the calculation
    print(f"Total area of the bounding square: {total_area}")
    print(f"Total number of samples: {num_samples}")
    print(f"Number of samples found in region R: {points_in_region}")
    print("\nThe area of R is estimated by the formula: (Points in Region / Total Samples) * Total Area")
    print(f"Area(R) = ({points_in_region} / {num_samples}) * {total_area}")
    print(f"The calculated area of R is: {estimated_area}")
    
    # Final answer formatted to two decimals
    print(f"\nFinal Answer (to two decimal places): {estimated_area:.2f}")

# Execute the calculation
calculate_area_of_R()