import math

def calculate_area_of_R():
    """
    Calculates the area of the region R defined by floor(|z|) = |floor(z)|
    within the square from 0 to 6+6i.
    """

    def G(t, R):
        """Helper function: Integral of sqrt(R^2 - t^2) from 0 to t."""
        if R <= 0:
            return 0.0
        # Clamp t to the domain of asin to handle floating point inaccuracies
        if abs(t) > R:
            t = math.copysign(R, t)
        
        R2 = R * R
        t2 = t * t
        
        val_sqrt = math.sqrt(max(0, R2 - t2))
        val_asin = math.asin(t / R)
        
        return (t * val_sqrt / 2) + (R2 * val_asin / 2)

    def calc_area_rect_disk(x, y, R):
        """
        Calculates the area of intersection of a rectangle [0,x]x[0,y]
        and a disk x^2+y^2<R^2 centered at the origin.
        This is done by integrating min(y, sqrt(R^2-t^2)) from t=0 to x.
        """
        if x <= 0 or y <= 0 or R <= 0:
            return 0.0
        
        # If the far corner of the rectangle is inside the circle, the area is the rectangle's area.
        if x * x + y * y <= R * R:
            return x * y

        R2 = R * R
        # tc is the x-coordinate where the circle intersects the line u=y
        tc = math.sqrt(max(0, R2 - y * y))

        if x <= tc:
            # In the integration range [0, x], y is always smaller than or equal to the circle curve
            return x * y
        else:
            # Integration range crosses tc, so the integral is split.
            # Area = integral from 0 to tc of y dt + integral from tc to x of sqrt(R^2-t^2) dt
            area = y * tc
            area += G(min(x, R), R) - G(tc, R)
            return area

    def area_square_disk_intersection(nx, ny, R):
        """
        Calculates the intersection area of a unit square [nx,nx+1]x[ny,ny+1] and a disk
        using the principle of inclusion-exclusion.
        """
        return (calc_area_rect_disk(nx + 1, ny + 1, R) -
                calc_area_rect_disk(nx, ny + 1, R) -
                calc_area_rect_disk(nx + 1, ny, R) +
                calc_area_rect_disk(nx, ny, R))

    total_area = 0.0
    
    # Find all integer pairs (nx, ny) in [0,5]x[0,5] where sqrt(nx^2+ny^2) is an integer C.
    valid_pairs = []
    for nx in range(6):
        for ny in range(6):
            C_squared = nx * nx + ny * ny
            C = math.sqrt(C_squared)
            # Check if C is an integer using a small tolerance for floating point comparison
            if abs(C - round(C)) < 1e-9:
                valid_pairs.append((nx, ny, int(round(C))))

    # Calculate area contribution from each valid unit square
    contributions = []
    print("The total area is the sum of contributions from each valid unit square [nx, nx+1] x [ny, ny+1].")
    print("The contribution from each square is the area of its intersection with the annulus C <= |z| < C+1.\n")
    
    print("Area =")
    for nx, ny, C in valid_pairs:
        # Area of square intersecting with outer disk (radius C+1)
        area_outer = area_square_disk_intersection(nx, ny, C + 1)
        # Area of square intersecting with inner disk (radius C)
        area_inner = area_square_disk_intersection(nx, ny, C)
        
        contribution = area_outer - area_inner
        total_area += contribution
        contributions.append(contribution)

    # Print the equation with all the numbers
    equation_parts = [f"{c:.4f}" for c in contributions]
    print(" + ".join(equation_parts))
    print(f"\nSumming these values, the total area is {total_area:.4f}.")
    print(f"The area of R, expressed to two decimals, is {round(total_area, 2)}.")

calculate_area_of_R()