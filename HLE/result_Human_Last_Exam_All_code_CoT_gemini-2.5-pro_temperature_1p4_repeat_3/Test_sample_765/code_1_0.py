import numpy as np
from scipy.integrate import dblquad

def calculate_tile_area(n, m):
    """
    Calculates the area of the valid region within a single unit square [n, n+1) x [m, m+1).
    """
    d_sq = n**2 + m**2
    d = np.sqrt(d_sq)
    
    # Return 0 if n^2 + m^2 is not a perfect square
    if d != np.floor(d):
        return 0.0
    
    K = int(d)
    # The condition on the region is K^2 <= x^2+y^2 < (K+1)^2
    # The lower bound is always satisfied for x>=n, y>=m.
    # We need to find the area of the square intersected with x^2+y^2 < (K+1)^2
    R_sq = (K + 1)**2
    
    # Define integration limits for the double integral
    # ∫[from x=n to n+1] ∫[from y=m to h(x)] 1 dy dx
    x_min, x_max = float(n), float(n + 1)
    
    # The upper limit for y is the lower of the square's top edge (m+1)
    # and the circle's boundary.
    y_lower_bound = lambda x: float(m)
    y_upper_bound = lambda x: min(m + 1.0, np.sqrt(R_sq - x**2)) if R_sq >= x**2 else float(m)
    
    # If the circle is entirely below the square, the area is 0.
    if R_sq <= m**2:
        return 0.0
    x_intersect_m = np.sqrt(R_sq - m**2)
    if x_intersect_m < n:
        return 0.0

    area, _ = dblquad(lambda y, x: 1.0, x_min, x_max, y_lower_bound, y_upper_bound)
    return area

if __name__ == '__main__':
    total_area = 0.0
    area_components = []
    
    # Find all pairs (n,m) where n^2+m^2 is a perfect square
    # and sort them for a structured output.
    valid_pairs = []
    for n in range(6):
        for m in range(6):
            d_sq = n**2 + m**2
            d = np.sqrt(d_sq)
            if d == np.floor(d):
                valid_pairs.append((n, m))
    
    # Calculate area for each valid pair
    for n, m in valid_pairs:
        area = calculate_tile_area(n, m)
        if area > 1e-6: # Filter out negligible areas
            total_area += area
            area_components.append(f"{area:.2f}")

    # Print the results as requested
    equation_str = " + ".join(area_components)
    print("The final area is a sum of contributions from different unit squares.")
    print("Final Equation:")
    print(f"{equation_str} = {total_area:.2f}")
    
    print(f"\nThe total area of R is {total_area:.2f}.")
