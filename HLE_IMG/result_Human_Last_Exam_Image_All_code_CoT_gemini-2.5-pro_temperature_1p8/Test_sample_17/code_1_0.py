import math

def main():
    """
    This program demonstrates that the theoretical stress at a perfectly sharp
    corner in a stressed elastic material is infinite.
    """

    # --- Assumptions for Demonstration ---
    # We assume some arbitrary but plausible values for the nominal stress and notch depth.
    # The final conclusion (infinite stress) does not depend on these specific values.
    sigma_nominal = 150.0  # Applied nominal stress in MPa (an example value for σ_y)
    notch_depth_a = 10.0   # Depth of the triangular wedge in mm

    print("This script calculates the theoretical stress at the tip of a sharp wedge (Point A).\n")
    print(f"Assumed Nominal Stress (σ_y) = {sigma_nominal} MPa")
    print(f"Assumed Wedge Depth (a) = {notch_depth_a} mm\n")

    print("The maximum stress is given by: σ_max = Kt * σ_nominal")
    print("The stress concentration factor (Kt) for a notch is approximated by: Kt = 1 + 2 * (a / ρ)")
    print("where 'ρ' is the radius of curvature at the notch tip.\n")
    print("Let's observe what happens as the tip gets sharper (ρ approaches 0):\n")

    # A list of progressively smaller tip radii to demonstrate the effect
    tip_radii_rho = [1.0, 0.1, 0.01, 0.001, 0.00001]

    for rho in tip_radii_rho:
        # Calculate Kt using the formula
        kt = 1 + 2 * (notch_depth_a / rho)
        # Calculate the maximum stress
        sigma_max = kt * sigma_nominal
        
        print(f"For a tip radius ρ = {rho:.5f} mm:")
        print(f"  Kt = 1 + 2 * ({notch_depth_a} / {rho:.5f}) = {kt:,.2f}")
        print(f"  σ_max = {kt:,.2f} * {sigma_nominal} = {sigma_max:,.2f} MPa")
        print("-" * 50)

    print("For a theoretically perfect sharp corner, the tip radius ρ is 0.")
    print(f"As ρ approaches 0, the term ({notch_depth_a} / ρ) approaches infinity.")
    print("Therefore, the stress concentration factor Kt also approaches infinity.")
    print("This means the theoretical stress at point A, σ_max, is infinite.")
    print("\nFinal Equation for the ideal case:")
    print("σ_max = (1 + 2 * (10.0 / 0)) * 150.0 = ∞")


if __name__ == "__main__":
    main()
