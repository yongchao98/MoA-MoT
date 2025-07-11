import math

def g(x, R):
    """
    Calculates the integral of sqrt(R^2 - t^2) from 0 to x.
    This represents the area of a sector and a triangle.
    """
    if R <= 0:
        return 0
    # Clamp x to the valid range [-R, R] for asin
    x_clamped = max(-R, min(R, x))
    return (x_clamped * math.sqrt(R**2 - x_clamped**2) / 2) + (R**2 * math.asin(x_clamped / R) / 2)

def area_circle_segment(R, x1, x2):
    """
    Calculates the area under the upper semi-circle of radius R from x1 to x2.
    The area is calculated as integral from x1 to x2 of sqrt(R^2 - x^2) dx.
    """
    if x1 >= x2 or R <= 0:
        return 0
    return g(x2, R) - g(x1, R)

def get_area_of_circle_in_square(na, nb, R):
    """
    Calculates the area of the circle x^2+y^2 < R^2 that lies within the
    unit square defined by [na, na+1] x [nb, nb+1].
    """
    # If the circle doesn't reach the square
    if R**2 <= na**2 + nb**2:
        return 0
    # If the square is completely inside the circle
    if R**2 >= (na + 1)**2 + (nb + 1)**2:
        return 1.0

    # Integration range for x
    x_start = na
    x_end = min(na + 1, R)

    # x-coordinates where the circle intersects the horizontal boundaries of the square
    x_intersect_top = -1
    if R > nb + 1:
        x_intersect_top = math.sqrt(R**2 - (nb + 1)**2)

    x_intersect_bottom = -1
    if R > nb:
        x_intersect_bottom = math.sqrt(R**2 - nb**2)

    # Collect all critical x-points for piecewise integration
    crit_x = sorted(list(set([x_start, x_end, x_intersect_top, x_intersect_bottom])))
    
    total_area = 0
    # Integrate over each segment defined by the critical points
    for i in range(len(crit_x) - 1):
        a = crit_x[i]
        b = crit_x[i+1]

        # Consider only the part of the interval within the square's x-range
        a = max(a, x_start)
        b = min(b, x_end)

        if a >= b:
            continue

        # Use a midpoint to determine the circle's position relative to the square
        mid_x = (a + b) / 2
        if R**2 < mid_x**2:
            continue
        y_circle = math.sqrt(R**2 - mid_x**2)

        if y_circle >= nb + 1:
            # Circle is above the square, integrand is (nb+1 - nb) = 1
            total_area += (b - a)
        elif y_circle > nb:
            # Circle crosses the square, integrand is sqrt(R^2-x^2) - nb
            total_area += area_circle_segment(R, a, b) - nb * (b - a)
        # else: circle is below the square, integrand is 0

    return total_area

def calculate_total_area():
    """
    Calculates the total area of the region R.
    """
    total_area = 0
    print("Calculating area contributions for each valid unit square [na, nb]:")
    print("-" * 60)
    print(f"{'Pair (na, nb)':<15}{'m':<5}{'Area Contribution':<20}")
    print("-" * 60)

    # Iterate over all unit squares in the 6x6 grid
    for na in range(6):
        for nb in range(6):
            s = na**2 + nb**2
            m_float = math.sqrt(s)
            
            # Check if na^2+nb^2 is a perfect square
            if m_float == math.floor(m_float):
                m = int(m_float)
                R1 = float(m)
                R2 = float(m + 1)
                
                # Area of the annulus in the square is Area(C(R2)) - Area(C(R1))
                area_R2 = get_area_of_circle_in_square(na, nb, R2)
                area_R1 = get_area_of_circle_in_square(na, nb, R1)
                
                contribution = area_R2 - area_R1
                
                if contribution > 1e-9: # Print only non-zero contributions
                    print(f"({na}, {nb}){'':<12}{m:<5}{contribution:<20.4f}")
                    total_area += contribution

    print("-" * 60)
    print(f"Total Area = {total_area}")
    print(f"The area of R is {total_area:.2f}")

if __name__ == '__main__':
    calculate_total_area()