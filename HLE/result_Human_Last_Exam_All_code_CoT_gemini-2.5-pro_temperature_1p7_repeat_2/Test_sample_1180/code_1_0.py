import math

def v(n):
    """
    Computes the 2-adic valuation of an integer n.
    Returns infinity for n=0.
    """
    if n == 0:
        return float('inf')
    if n % 2 != 0:
        return 0
    return (n & -n).bit_length() - 1

def main():
    """
    Calculates the thickness of the double point for the given curve.
    """
    # The coefficients of the polynomial q(u) = 64u^5 + 160u^4 + 176u^3 + 104u^2 + 32u + 5
    coeffs = [5, 32, 104, 176, 160, 64]
    
    # Calculate the 2-adic valuations of the coefficients
    valuations = [v(c) for c in coeffs]
    
    # The points for the Newton polygon are (i, v(c_i))
    points = list(enumerate(valuations))
    
    print("The polynomial in u is: q(u) = 64*u^5 + 160*u^4 + 176*u^3 + 104*u^2 + 32*u + 5")
    print("The points for the Newton polygon of q(u) are (i, v(c_i)):")
    for i, p in enumerate(points):
        print(f"  For u^{i}: point is {p}")
        
    # The vertices of the lower convex hull are (0,0), (2,3), and (5,6).
    # We are interested in the first segment, corresponding to the two close roots.
    p1 = points[0]  # (0, 0)
    p2 = points[2]  # (2, 3)
    
    # Calculate the slope of the first segment
    slope_num = p2[1] - p1[1]
    slope_den = p2[0] - p1[0]
    
    thickness = slope_num / slope_den
    
    print("\nThe lower convex hull of these points has vertices at (0,0), (2,3), and (5,6).")
    print("The double point corresponds to the two roots that are close to each other.")
    print("The thickness of this double point is the slope of the first segment of the Newton polygon.")
    print(f"Slope = (v(c_2) - v(c_0)) / (2 - 0) = ({p2[1]} - {p1[1]}) / ({p2[0]} - {p1[0]}) = {slope_num}/{slope_den}")
    print(f"The thickness is: {thickness}")

if __name__ == "__main__":
    main()
