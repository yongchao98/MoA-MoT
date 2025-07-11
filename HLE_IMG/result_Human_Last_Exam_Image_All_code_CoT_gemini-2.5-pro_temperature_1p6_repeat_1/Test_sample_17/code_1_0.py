import math

# --- User-defined variables for demonstration ---
# Let's assume some values for the problem.

# Nominal uniaxial stress applied to the plate (e.g., in MPa)
sigma_y = 100

# Depth of the triangular wedge/notch (e.g., in mm)
d = 5.0

# Tip radius at point A (e.g., in mm).
# We'll use a very small, non-zero value to demonstrate the effect.
# In theory, for a perfect corner, rho = 0.
rho = 0.0001

# --- Calculation ---
# For a V-notch, a common approximation for the stress concentration factor (Kt) is:
# Kt = 1 + 2 * sqrt(d / rho)
# This shows that as rho -> 0, Kt -> infinity.

if rho > 0:
    Kt = 1 + 2 * math.sqrt(d / rho)
    sigma_max = Kt * sigma_y

    # --- Output the results ---
    print(f"Demonstration of Stress Concentration at a Notch Tip")
    print("-" * 50)
    print(f"Given values:")
    print(f"  Nominal Stress (σ_y): {sigma_y} MPa")
    print(f"  Notch Depth (d): {d} mm")
    print(f"  Tip Radius (ρ): {rho} mm (a very small, non-zero value)")
    print("-" * 50)
    print("Step 1: Calculate the Stress Concentration Factor (Kt)")
    print(f"  Kt = 1 + 2 * sqrt(d / ρ)")
    print(f"  Kt = 1 + 2 * sqrt({d} / {rho})")
    print(f"  Kt = 1 + 2 * sqrt({d/rho})")
    print(f"  Kt = 1 + 2 * {math.sqrt(d/rho):.2f}")
    print(f"  Kt = {Kt:.2f}")
    print("")
    print("Step 2: Calculate the Maximum Stress at the Tip (σ_max)")
    print(f"  σ_max = Kt * σ_y")
    print(f"  σ_max = {Kt:.2f} * {sigma_y}")
    print(f"  σ_max = {sigma_max:.2f} MPa")
    print("-" * 50)
    print("Conclusion: As the tip radius 'ρ' approaches zero, 'Kt' and therefore 'σ_max' approach infinity.")
    print("The theoretical stress at a perfectly sharp tip is infinite.")
else:
    print("The theoretical stress at a perfectly sharp tip (rho=0) is infinite.")
