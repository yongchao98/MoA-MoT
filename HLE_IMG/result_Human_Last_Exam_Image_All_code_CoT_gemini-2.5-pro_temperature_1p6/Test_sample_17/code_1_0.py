import math

def calculate_theoretical_stress_at_wedge_tip():
    """
    Calculates and explains the theoretical stress at the tip of a sharp wedge
    in a plate under uniaxial tension.
    """
    # The geometry shown, with sharp re-entrant corners (the wedge tips),
    # creates a stress concentration. According to the theory of linear elasticity,
    # for a perfectly sharp corner, the stress at the very tip is singular.

    # The stress field near the tip of a sharp notch (point A) can be approximated
    # by the equations for a crack tip. The stress 'σ' at a small distance 'r'
    # from the tip is given by:
    # σ(r) = K_I / sqrt(2 * π * r)
    # Where K_I is the stress intensity factor, which depends on the applied
    # stress (σ_y) and the geometry.

    print("The theoretical stress 'σ' at a small distance 'r' from the sharp tip (Point A) is given by the equation:")
    print("σ(r) = K_I / sqrt(2 * π * r)")
    print("-" * 30)

    # To find the stress *at* the tip, we must evaluate the limit as r -> 0.
    print("The stress *at* the tip itself is the limit of this equation as r approaches 0.")

    # Displaying the final equation leading to the result.
    # The numbers in the equation are 2 and π.
    pi_val = math.pi
    print("\nFinal Equation to be Evaluated:")
    print("σ_tip = lim (r→0) [ K_I / sqrt(2 * π * r) ]")
    print(f"σ_tip = lim (r→0) [ K_I / sqrt({2} * {pi_val:.4f} * r) ]")
    
    # As r approaches 0, the denominator sqrt(2 * pi * r) also approaches 0.
    # A finite value (K_I) divided by a value approaching zero results in infinity.
    result = "Infinity"
    
    print("\nResult:")
    print(f"As the distance 'r' approaches 0, the denominator approaches 0, causing the expression to go to infinity.")
    print(f"Therefore, the theoretical stress at the tip is: {result}")
    print("\n(Note: This is a theoretical result. In real materials, stress is limited by the material's yield strength, leading to local plastic deformation instead of infinite stress.)")

calculate_theoretical_stress_at_wedge_tip()
<<<Infinity>>>