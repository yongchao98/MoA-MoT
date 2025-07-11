import math

def calculate_stress_at_notch_tip(sigma_nominal, notch_depth_a, tip_radius_rho):
    """
    Calculates the theoretical maximum stress at the tip of a notch
    using the formula for an elliptical notch to demonstrate stress singularity.

    Args:
        sigma_nominal (float): The nominal applied stress (σ_y).
        notch_depth_a (float): The depth of the notch.
        tip_radius_rho (float): The radius of curvature at the notch tip.

    Returns:
        tuple(float, float): A tuple containing the maximum theoretical stress
                             at the notch tip and the stress concentration factor.
    """
    if tip_radius_rho <= 0:
        return float('inf'), float('inf')

    # Stress concentration factor Kt = 1 + 2 * sqrt(a / rho)
    kt = 1 + 2 * math.sqrt(notch_depth_a / tip_radius_rho)

    # Maximum stress sigma_max = Kt * sigma_nominal
    sigma_max = kt * sigma_nominal
    return sigma_max, kt

# --- Main Program ---

# Define placeholder variables for the purpose of demonstration.
# The applied stress σ_y is the nominal stress in the far-field.
sigma_y = 100.0  # Assumed nominal applied stress in MPa
a = 2.0          # Assumed notch depth in mm

print("This program demonstrates the concept of a stress singularity at a sharp notch (Point A).")
print("According to linear elastic theory, the stress at a perfectly sharp corner (tip radius ρ -> 0) is infinite.")
print("We can demonstrate this using the stress concentration formula for an elliptical notch:")
print("σ_max = σ_nominal * (1 + 2 * sqrt(a/ρ))")
print("-" * 70)
print(f"Assuming a Nominal Stress σ_y = {sigma_y} MPa and a Notch Depth a = {a} mm.")
print("Let's observe σ_max as the tip radius ρ gets closer to zero:")
print("-" * 70)

# Demonstrate with progressively smaller (sharper) tip radii
tip_radii = [1e-1, 1e-3, 1e-6, 1e-9, 1e-12]

for rho in tip_radii:
    sigma_max, kt = calculate_stress_at_notch_tip(sigma_y, a, rho)
    print(f"For a tip radius ρ = {rho:.1e} mm:")
    
    # Printing the calculation for the stress concentration factor Kt
    kt_calculation_str = f"Kt = 1 + 2 * sqrt({a} / {rho:.1e})"
    print(f"  {kt_calculation_str} = {kt:,.2f}")
    
    # Printing the final equation for sigma_max with all its numbers
    sigma_max_equation_str = f"σ_max = {sigma_y} * {kt:,.2f}"
    print(f"  {sigma_max_equation_str} = {sigma_max:,.2f} MPa")
    print()

print("Conclusion:")
print("As the notch tip becomes perfectly sharp (ρ -> 0), the stress concentration")
print("factor (Kt) and the maximum stress (σ_max) grow without bound.")
print("\nTherefore, the THEORETICAL stress value at the tip of the wedge (Point A) is infinite.")
