import math
from fractions import Fraction

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem near the ground
    using the mirror image method.
    """
    # 1. Define geometric parameters. We can set c=1 as it will cancel out.
    c = 1.0
    h = c / 2
    s = c / 2

    # 2. Define vortex and control point locations.
    # Vortex is at c/4, Control Point is at 3c/4.
    # Aerofoil 1 (front): Placed with its leading edge at x=0.
    x_v1 = c / 4
    y_v1 = h
    x_c1 = 3 * c / 4
    y_c1 = h

    # Aerofoil 2 (rear): Placed with its leading edge at x = c + s.
    x_v2 = (c + s) + (c / 4)
    y_v2 = h
    x_c2 = (c + s) + (3 * c / 4)
    y_c2 = h
    
    # 3. Define a function to calculate the vertical velocity 'w' induced by a vortex 'Gamma'.
    # The coefficient is w/Gamma.
    def get_w_coeff(x_v, y_v, x_c, y_c):
        """Calculates the coefficient for induced vertical velocity."""
        dx = x_c - x_v
        dy = y_c - y_v
        r_squared = dx**2 + dy**2
        if r_squared < 1e-9:  # Avoid division by zero
            return 0
        return (1 / (2 * math.pi)) * (dx / r_squared)

    # 4. Calculate the influence coefficients C_ij for the system w_i = sum(C_ij * Gamma_j).
    # Note: The influence of a vortex on its own aerofoil is captured by the main lift equation.
    # The coefficients C_ij account for interference from *other* vortices (including images).

    # C11: Influence on P1 from image of Gamma1 (-Gamma1)
    C11 = -1 * get_w_coeff(x_v1, -y_v1, x_c1, y_c1)
    # C12: Influence on P1 from Gamma2 and its image (-Gamma2)
    C12 = get_w_coeff(x_v2, y_v2, x_c1, y_c1) - get_w_coeff(x_v2, -y_v2, x_c1, y_c1)

    # C21: Influence on P2 from Gamma1 and its image (-Gamma1)
    C21 = get_w_coeff(x_v1, y_v1, x_c2, y_c2) - get_w_coeff(x_v1, -y_v1, x_c2, y_c2)
    # C22: Influence on P2 from image of Gamma2 (-Gamma2)
    C22 = -1 * get_w_coeff(x_v2, -y_v2, x_c2, y_c2)

    # 5. Form the system of equations.
    # The governing equation for each aerofoil is: Gamma_i = pi*c*U*alpha + pi*c*w_i
    # Let K = pi*c and rhs = K*U*alpha (which is the same for both).
    #
    # (1 - K*C11) * Gamma1 - (K*C12) * Gamma2 = rhs
    # -(K*C21) * Gamma1 + (1 - K*C22) * Gamma2 = rhs
    
    K = math.pi * c
    # Matrix A coefficients for A * [Gamma1, Gamma2]^T = B
    A11 = 1 - K * C11
    A12 = -K * C12
    A21 = -K * C21
    A22 = 1 - K * C22

    # Since the RHS is the same, we can equate the two left-hand sides:
    # A11*Gamma1 + A12*Gamma2 = A21*Gamma1 + A22*Gamma2
    # (A11 - A21)*Gamma1 = (A22 - A12)*Gamma2
    coeff_G1 = A11 - A21
    coeff_G2 = A22 - A12
    
    # 6. Present the final equation and solve for the ratio L1/L2.
    # The relation coeff_G1 * Gamma1 = coeff_G2 * Gamma2 is equivalent to
    # coeff_G1 * L1 = coeff_G2 * L2. We want to find integer coefficients.
    f_G1 = Fraction(coeff_G1).limit_denominator(100)
    f_G2 = Fraction(coeff_G2).limit_denominator(100)

    # Simplify to get X * L1 = Y * L2 where X, Y are coprime integers.
    common_divisor = math.gcd(f_G1.numerator * f_G2.denominator, f_G2.numerator * f_G1.denominator)
    X = (f_G1.numerator * f_G2.denominator) // common_divisor
    Y = (f_G2.numerator * f_G1.denominator) // common_divisor

    print("The final relationship between the lift on the two aerofoils (L1 and L2) is:")
    print(f"{X} * L1 = {Y} * L2")

    lift_ratio = Y / X
    print("\nThe lift ratio L1/L2 is:")
    print(f"{Y}/{X}")
    print(f"= {lift_ratio}")
    return lift_ratio

# Run the calculation and store the result
final_ratio = calculate_lift_ratio()
print(f"\n<<<solution\n{final_ratio}\n>>>")