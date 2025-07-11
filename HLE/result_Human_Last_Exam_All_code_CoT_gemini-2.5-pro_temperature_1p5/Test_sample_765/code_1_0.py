import math

def calculate_area():
    """
    Calculates the area of the region R using numerical integration.
    """
    
    # The overall domain is a 6x6 square
    domain_min, domain_max = 0, 6
    
    # Find all integer pairs (nx, ny) in the range [0, 5] x [0, 5]
    # for which sqrt(nx^2 + ny^2) is an integer.
    valid_squares = []
    for nx in range(domain_max):
        for ny in range(domain_max):
            k_float = math.sqrt(nx*nx + ny*ny)
            # Check if k_float is an integer (within a small tolerance for floating point)
            if abs(k_float - round(k_float)) < 1e-9:
                valid_squares.append((nx, ny))

    # Step size for the grid. A smaller delta increases accuracy.
    delta = 0.002
    cell_area = delta * delta
    
    total_area = 0.0
    area_terms = []

    # Iterate over each valid unit square and calculate its contribution to the area
    for nx, ny in valid_squares:
        k = int(round(math.sqrt(nx*nx + ny*ny)))
        square_area = 0.0
        
        # Define the boundaries of the current unit square
        x_start, x_end = nx, nx + 1
        y_start, y_end = ny, ny + 1
        
        # Iterate over the grid points within the unit square
        num_steps = int(1.0 / delta)
        for i in range(num_steps):
            x = x_start + (i + 0.5) * delta
            for j in range(num_steps):
                y = y_start + (j + 0.5) * delta
                
                # Check if the point satisfies the condition for the region R
                # The condition for this square is k <= |z| < k+1
                mod_z_sq = x*x + y*y
                if k*k <= mod_z_sq < (k+1)*(k+1):
                    square_area += cell_area
                    
        area_terms.append(square_area)
        total_area += square_area

    # Output the results as requested
    print("The total area is the sum of the areas from each valid unit square:")
    equation_parts = [f"{term:.2f}" for term in area_terms]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_area:.2f}")

    print(f"\nThe area of the region R is approximately: {total_area:.2f}")

if __name__ == '__main__':
    calculate_area()
