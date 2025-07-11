import numpy as np

def polygon_area(x, y):
    """
    Calculates the area of a polygon using the Shoelace Formula.
    The vertices (x, y) are given as two lists.
    """
    # The formula requires looping back to the first vertex
    x_shifted = np.roll(x, -1)
    y_shifted = np.roll(y, -1)
    
    # Shoelace formula components
    term1 = np.sum(x * y_shifted)
    term2 = np.sum(y * x_shifted)
    
    area = 0.5 * np.abs(term1 - term2)
    return area

def solve_minimal_area():
    """
    Solves for the minimal area of a convex domain intersecting all lines
    px+qy=1 for coprime integers p and q.
    
    The minimal domain is a specific hexagon. This function calculates its area.
    """
    # Vertices of the hexagon
    # H = conv{(1,0), (0,1), (-1,1), (-1,0), (0,-1), (1,-1)}
    x_coords = np.array([1, 0, -1, -1, 0, 1])
    y_coords = np.array([0, 1,  1,  0, -1,-1])
    
    # Calculate the area
    area = polygon_area(x_coords, y_coords)
    
    # The result is used in an equation to show the calculation.
    # The equation is the Shoelace formula for this hexagon.
    # A = 0.5 * |(x0y1 + x1y2 + ... ) - (y0x1 + y1x2 + ...)|
    term1_vals = x_coords * np.roll(y_coords, -1)
    term2_vals = y_coords * np.roll(x_coords, -1)
    
    print("The minimal area is the area of a hexagon with vertices at (1,0), (0,1), (-1,1), (-1,0), (0,-1), and (1,-1).")
    print("We can calculate this area using the Shoelace formula:")
    print("Area = 0.5 * |(x0*y1 + x1*y2 + ...) - (y0*x1 + y1*x2 + ...)|")
    
    term1_str = " + ".join([f"({x}*{y})" for x, y in zip(x_coords, np.roll(y_coords, -1))])
    term2_str = " + ".join([f"({y}*{x})" for y, x in zip(y_coords, np.roll(x_coords, -1))])
    
    print(f"Area = 0.5 * |({term1_str}) - ({term2_str})|")
    
    sum1 = np.sum(term1_vals)
    sum2 = np.sum(term2_vals)
    
    print(f"Area = 0.5 * |({sum1}) - ({sum2})|")
    print(f"Area = 0.5 * |{sum1 - sum2}|")
    print(f"Area = {0.5 * abs(sum1 - sum2)}")
    
    print("\nThe minimal area is:")
    print(area)

solve_minimal_area()