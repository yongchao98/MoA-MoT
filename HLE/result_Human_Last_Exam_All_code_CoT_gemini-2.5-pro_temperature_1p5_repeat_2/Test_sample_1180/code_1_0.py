import math

def get_valuation(n):
    """Computes the 2-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    if n % 2 != 0:
        return 0
    return (n & -n).bit_length() - 1

def main():
    """
    Calculates the thickness of the double point of the stable reduction of the curve
    z^2 = 2*x^5 + 2*x^3 + 1.
    """
    # The curve is z^2 = 2*x^5 + 2*x^3 + 1.
    # The coefficients of the polynomial f(x) are a_5=2, a_3=2, a_0=1.
    coeffs = {5: 2, 3: 2, 0: 1}

    # The points for the Newton polygon are (i, v(a_i)).
    # v(2)=1, v(1)=0.
    points = sorted([(i, get_valuation(coeffs[i])) for i in coeffs])

    # The points are [(0, 0), (3, 1), (5, 1)].
    # We compute the lower convex hull and its slopes.
    # The vertices of the lower convex hull are (0,0), (3,1), and (5,1).
    p1 = points[0]  # (0, 0)
    p2 = points[1]  # (3, 1)
    p3 = points[2]  # (5, 1)

    # Slope between p2 and p1
    dy = p2[1] - p1[1]
    dx = p2[0] - p1[0]
    
    # The denominator of this slope is interpreted as the thickness.
    # The fraction is dy/dx. We need its denominator.
    common_divisor = math.gcd(dy, dx)
    numerator = dy // common_divisor
    denominator = dx // common_divisor
    
    print("The Newton Polygon is determined by the points (i, v(a_i)).")
    print(f"The relevant points for the slope calculation are ({p1[0]}, {p1[1]}) and ({p2[0]}, {p2[1]}).")
    print("The slope is calculated as (y2 - y1) / (x2 - x1).")
    print(f"y2 - y1 = {p2[1]} - {p1[1]} = {dy}")
    print(f"x2 - x1 = {p2[0]} - {p1[0]} = {dx}")
    print(f"The resulting slope is {numerator}/{denominator}.")
    print("The thickness is the denominator of this fractional slope.")
    print(f"The thickness is {denominator}.")

main()