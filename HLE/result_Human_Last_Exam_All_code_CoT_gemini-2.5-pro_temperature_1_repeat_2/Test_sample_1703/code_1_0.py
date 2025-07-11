import math

def solve_triangle_area():
    """
    This function explains the derivation and prints the final formula for the area
    of the triangle T(t) as a function of time t.
    """
    
    # Define problem constants
    R = 10.0  # Radius of the circle and side length of the hexagon
    v = 1.0   # Speed of the triangle's vertices

    print("--- Derivation of the Area of Triangle T(t) ---")
    print("Step 1: Simplify the problem.")
    print("The rotation of the hexagon does not change the triangle's area. We can analyze it in a static frame.")
    
    print("\nStep 2: Use symmetry.")
    print("The setup is symmetric, so the triangle T(t) is always equilateral. We only need its side length, S(t).")
    print("The area will be Area = (sqrt(3)/4) * S(t)^2.")

    print("\nStep 3: Use coordinate geometry to find S(t)^2.")
    print(f"The hexagon has side length s = R = {R}.")
    s = R
    
    print("Let the vertices of T(t) be P(t) on side AB, Q(t) on side CD, and R(t) on side EF.")
    print("We will find the coordinates of P(t) and Q(t) to calculate the distance between them.")
    
    # In a static frame with center (0,0), the coordinates of P(t) and Q(t) can be derived.
    # P(t) moves on side AB from the midpoint. Its coordinates are:
    # Px(t) = (3*s - 2*t)/4 = (15 - t)/2
    # Py(t) = (s + 2*t)*sqrt(3)/4 = (5 + t)*sqrt(3)/2
    #
    # Q(t) moves on side CD from the midpoint. Its coordinates are:
    # Qx(t) = (-3*s - 2*t)/4 = (-15 - t)/2
    # Qy(t) = (s - 2*t)*sqrt(3)/4 = (5 - t)*sqrt(3)/2
    
    # The squared distance S(t)^2 = (Px - Qx)^2 + (Py - Qy)^2
    # Px - Qx = (15 - t)/2 - (-15 - t)/2 = 30 / 2 = 15
    # Py - Qy = ((5 + t) - (5 - t))*sqrt(3)/2 = 2*t*sqrt(3)/2 = t*sqrt(3)
    # S(t)^2 = 15^2 + (t*sqrt(3))^2
    
    side_len_sq_const = (3 * s / 2)**2
    side_len_sq_t_coeff = 3

    print(f"\nThe squared side length is calculated as S(t)^2 = {side_len_sq_const:.0f} + {side_len_sq_t_coeff}*t^2.")

    print("\nStep 4: Calculate the area.")
    print(f"Substituting S(t)^2 into the area formula for an equilateral triangle:")
    print(f"Area(t) = (sqrt(3)/4) * ({side_len_sq_const:.0f} + {side_len_sq_t_coeff}*t^2)")

    # The simplified form is obtained by factoring out 3.
    # Area(t) = (3 * sqrt(3) / 4) * (75 + t^2)
    num1 = 3
    num2_in_sqrt = 3
    den = 4
    num3 = int(side_len_sq_const / num1)

    print("\n--- Final Formula ---")
    print("The area of the triangle T(t) as a function of time t is:")
    print(f"Area(t) = ({num1} * sqrt({num2_in_sqrt}) / {den}) * ({num3} + t^2)")

if __name__ == '__main__':
    solve_triangle_area()