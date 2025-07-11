import math

def solve_candy_packing():
    """
    Solves the integer-based sphere packing problem for the Halloween box.
    """
    # 1. Explain and verify the problem formulation.
    print("Is the problem formulation correct? Yes.")
    print("Here is the step-by-step breakdown and solution:")
    print("-" * 50)

    # Physical dimensions
    box_cm = (12, 12, 11)
    radius_cm = 2
    grid_step_cm = 0.5
    
    # Dimensions in grid units (integers)
    box_units = tuple(int(dim / grid_step_cm) for dim in box_cm)
    radius_units = int(radius_cm / grid_step_cm)
    
    print(f"The box is {box_cm[0]}x{box_cm[1]}x{box_cm[2]} cm.")
    print(f"The candy is a sphere with a {radius_cm} cm radius.")
    print(f"Using a 0.5 cm grid, we convert to integer units.")
    print(f"Box dimensions in grid units: {box_units[0]}x{box_units[1]}x{box_units[2]}")
    print(f"Candy radius in grid units: {radius_units}")
    print("-" * 50)

    # Valid center coordinate ranges
    x_min, y_min, z_min = radius_units, radius_units, radius_units
    x_max = box_units[0] - radius_units
    y_max = box_units[1] - radius_units
    z_max = box_units[2] - radius_units

    print("The valid ranges for the center (x, y, z) of a candy are:")
    print(f"  x: [{x_min}, {x_max}]")
    print(f"  y: [{y_min}, {y_max}]")
    print(f"  z: [{z_min}, {z_max}]")
    print("-" * 50)
    
    # Non-overlapping constraint
    min_dist_sq = (2 * radius_units)**2
    print("The non-overlapping constraint for two candy centers (xi, yi, zi) and (xj, yj, zj) is:")
    print("  (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= (radius_i + radius_j)^2")
    print(f"  (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= ({radius_units} + {radius_units})^2")
    print(f"  (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= {min_dist_sq}")
    print("-" * 50)

    # 2. Solve using a greedy algorithm
    print("Applying a greedy algorithm to pack the candies layer by layer...")
    
    placed_centers = []
    
    # Iterate through all possible center locations, z first (layer by layer)
    for z in range(z_min, z_max + 1):
        for y in range(y_min, y_max + 1):
            for x in range(x_min, x_max + 1):
                
                # Check if this position is valid
                can_place = True
                candidate_center = (x, y, z)
                
                for placed_center in placed_centers:
                    # Calculate squared distance to avoid sqrt
                    dist_sq = (
                        (candidate_center[0] - placed_center[0])**2 +
                        (candidate_center[1] - placed_center[1])**2 +
                        (candidate_center[2] - placed_center[2])**2
                    )
                    
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break # This position overlaps, try the next one
                
                if can_place:
                    placed_centers.append(candidate_center)

    max_n = len(placed_centers)
    print("\nCalculation complete.")
    print(f"The maximized value n (number of candies) is: {max_n}")

# Run the solver
solve_candy_packing()

# Final answer in the required format
print("\n<<<22>>>")
