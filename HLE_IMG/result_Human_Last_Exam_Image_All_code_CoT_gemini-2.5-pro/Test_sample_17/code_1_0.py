import math

def calculate_stress_at_notch_tip():
    """
    This function demonstrates the concept of stress singularity at a sharp corner.

    In theory, the stress at a perfectly sharp corner (like point A) in an elastic
    material is infinite. This is because the radius of curvature at the tip is
    zero, leading to an infinite stress concentration factor.

    We can demonstrate this by calculating the maximum stress for progressively
    smaller tip radii using a common formula for stress concentration at a notch:
    σ_max = σ_nom * (1 + 2 * sqrt(a / ρ))
    where:
    - σ_max is the maximum stress at the tip.
    - σ_nom is the nominal applied stress.
    - 'a' is the depth of the notch (or wedge).
    - 'ρ' (rho) is the radius of curvature at the notch tip.
    """

    # --- Assumptions for Demonstration ---
    # Let's assume a nominal applied stress (σ_y in the diagram).
    sigma_nom = 100.0  # in Megapascals (MPa)

    # Let's assume a depth 'a' for the triangular wedge.
    a = 2.0  # in millimeters (mm)

    # We will now observe the effect of making the tip sharper by decreasing the radius ρ.
    # Note: ρ is the radius of the very tip of the point, not the overall wedge geometry.
    rho_values = [0.1, 0.01, 0.001, 0.0001, 1e-6, 1e-9]

    print(f"Demonstrating Stress Concentration with:")
    print(f"Nominal Stress (σ_nom) = {sigma_nom} MPa")
    print(f"Notch Depth (a) = {a} mm\n")
    print("-" * 60)
    print(f"{'Tip Radius ρ (mm)':<20} | {'Stress Concentration Kt':<25} | {'Max Stress σ_max (MPa)':<25}")
    print("-" * 60)

    for rho in rho_values:
        # Calculate the Stress Concentration Factor (Kt)
        # Equation: Kt = 1 + 2 * sqrt(a / ρ)
        kt = 1 + 2 * math.sqrt(a / rho)

        # Calculate the maximum stress (σ_max)
        # Equation: σ_max = σ_nom * Kt
        sigma_max = sigma_nom * kt

        # Print the results showing the equation with numbers
        print(f"{rho:<20.1e} | 1 + 2*sqrt({a}/{rho:.1e}) = {kt:<10.2e} | {sigma_nom} * {kt:.2e} = {sigma_max:<25.2e}")

    print("-" * 60)
    print("\nConclusion:")
    print("As the tip radius 'ρ' approaches 0, the stress concentration factor 'Kt' and the")
    print("maximum stress 'σ_max' approach infinity.")
    print("Therefore, the theoretical stress at the perfectly sharp tip (point A) is infinite.")

# Run the calculation and print the explanation.
calculate_stress_at_notch_tip()