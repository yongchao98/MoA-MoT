import math

def demonstrate_stress_singularity():
    """
    This function demonstrates the concept of a stress singularity at a sharp corner.
    """

    # Assume a nominal (far-field) stress applied to the plate.
    # The unit can be any stress unit, e.g., MPa (Megapascals).
    sigma_nominal = 100.0

    # Assume a characteristic dimension for the notch, its depth.
    # The unit must be consistent with the radius unit, e.g., mm.
    notch_depth = 2.0

    print("--- Stress Concentration at a Sharp Corner ---")
    print("The theoretical stress at a perfectly sharp corner (radius = 0) is infinite.")
    print("This is because the Stress Concentration Factor (Kt) approaches infinity as the tip radius approaches zero.")
    print("\nLet's demonstrate this with an approximate formula for Kt:")
    print("  σ_max = Kt * σ_nominal")
    print("  Kt ≈ 1 + 2 * sqrt(notch_depth / tip_radius)")
    print("--------------------------------------------------")
    print(f"Assuming Nominal Stress = {sigma_nominal} MPa, Notch Depth = {notch_depth} mm\n")

    # A list of decreasing tip radii to simulate approaching a sharp corner
    tip_radii = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.000001]

    for radius in tip_radii:
        # Calculate the Stress Concentration Factor (Kt)
        kt = 1 + 2 * math.sqrt(notch_depth / radius)
        # Calculate the maximum stress at the tip
        sigma_max = kt * sigma_nominal

        print(f"For a tip radius of {radius:.6f} mm:")
        # Displaying the numbers in the final equation as requested
        print(f"  Final Equation: σ_max = ({1:.1f} + {2:.1f} * sqrt({notch_depth:.1f} / {radius:.6f})) * {sigma_nominal:.1f}")
        print(f"  Result:         σ_max = {kt:.2f} * {sigma_nominal:.1f} = {sigma_max:,.2f} MPa\n")

    print("Conclusion: As the tip radius approaches zero, the calculated maximum stress approaches infinity.")
    print("Therefore, the theoretical stress at point A is infinite.")

if __name__ == "__main__":
    demonstrate_stress_singularity()