import numpy as np

def demonstrate_stress_singularity(nominal_stress, notch_depth):
    """
    Demonstrates that theoretical stress at a sharp corner is infinite.

    Args:
        nominal_stress (float): The far-field stress applied to the plate (sigma_y).
        notch_depth (float): The depth of the triangular notch (a).
    """

    print("This program demonstrates the concept of stress concentration at a sharp corner.")
    print(f"We use a nominal stress (sigma_y) of {nominal_stress} MPa and a notch depth (a) of {notch_depth} mm.")
    print("The maximum stress (sigma_max) is given by: sigma_max = K_t * sigma_y")
    print("Where the stress concentration factor (K_t) for a notch is approximately: K_t = 1 + 2 * sqrt(a / rho)")
    print("'rho' is the radius of the tip. For a perfectly sharp corner, rho is theoretically 0.\n")
    print("Let's observe the stress as 'rho' approaches 0:\n")

    # A list of progressively smaller tip radii
    tip_radii = [0.1, 0.01, 0.001, 0.0001]

    for rho in tip_radii:
        # Calculate the stress concentration factor
        Kt = 1 + 2 * np.sqrt(notch_depth / rho)
        # Calculate the maximum stress
        sigma_max = Kt * nominal_stress

        print(f"For a tip radius rho = {rho} mm:")
        # The final equation is sigma_max = nominal_stress * (1 + 2 * sqrt(notch_depth / rho))
        print(f"  sigma_max = {nominal_stress} * (1 + 2 * sqrt({notch_depth} / {rho}))")
        print(f"  sigma_max = {nominal_stress} * ({Kt:.2f}) = {sigma_max:.2f} MPa")
        print("-" * 30)

    print("As the tip radius 'rho' approaches 0, the stress concentration factor K_t and the maximum stress 'sigma_max' approach infinity.")
    print("\nTherefore, the theoretical stress at the tip of the wedge (point A) is infinite.")

# --- Parameters for Demonstration ---
# Let's assume some reasonable values for nominal stress and notch depth.
sigma_y_nominal = 150.0  # MPa
notch_depth_a = 1.0      # mm

# Run the demonstration
demonstrate_stress_singularity(sigma_y_nominal, notch_depth_a)