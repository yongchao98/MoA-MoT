import numpy as np
from scipy.integrate import quad

def integrand(x, y1, y2, R):
    """
    Calculates the length of the vertical slice of the intersection
    of a disk x^2+y^2 < R^2 and a rectangle [x1, x2] x [y1, y2] at a given x.
    This function is designed to be integrated over x from x1 to x2.
    """
    # If x is outside the circle's radius, the slice length is 0.
    if x >= R:
        return 0.0
    
    # y-coordinate on the upper semi-circle
    y_circle = np.sqrt(R**2 - x**2)
    
    # The height of the intersection slice is determined by the rectangle's
    # y-boundaries and the circle's boundary.
    height = max(0, min(y2, y_circle) - y1)
    return height

def calculate_intersect_area(na, nb, R):
    """
    Calculates the area of intersection between the unit square
    [na, na+1] x [nb, nb+1] and the disk x^2 + y^2 < R^2.
    """
    x1, x2 = na, na + 1
    y1, y2 = nb, nb + 1
    
    # Integrate the slice height function over the width of the square.
    area, error = quad(integrand, x1, x2, args=(y1, y2, R))
    return area

def solve_area():
    """
    Calculates the total area of the region R.
    """
    total_area = 0
    area_terms = []

    # Iterate over all unit squares within the 6x6 grid.
    for na in range(6):
        for nb in range(6):
            k_float = np.sqrt(na**2 + nb**2)
            
            # Check if k = sqrt(na^2 + nb^2) is an integer.
            if np.isclose(k_float, round(k_float)):
                k = int(round(k_float))
                
                # For this square, the region is bounded by a circle of radius R = k + 1.
                R = k + 1
                
                # Calculate the area contribution from this square.
                area_sq = calculate_intersect_area(na, nb, R)
                total_area += area_sq
                area_terms.append(f"{area_sq:.2f}")

    # Print the equation representing the sum of the areas.
    print(f"The total area is the sum of the areas from each qualifying unit square:")
    print("Area = " + " + ".join(area_terms))
    print(f"Total Area = {total_area:.2f}")

solve_area()