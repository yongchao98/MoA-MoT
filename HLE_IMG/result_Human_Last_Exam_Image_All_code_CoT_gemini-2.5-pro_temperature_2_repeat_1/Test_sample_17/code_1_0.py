import sympy as sp
import math

def solve_stress_concentration():
    """
    Calculates and explains the theoretical stress at the tip of a sharp wedge.
    """

    # 1. Define symbolic variables
    sigma_nom = sp.Symbol('σ_nominal', positive=True) # Nominal applied stress
    d = sp.Symbol('d', positive=True)              # Depth of the notch/wedge
    rho = sp.Symbol('ρ', positive=True)             # Radius of curvature at the tip of the notch

    # 2. Explain the theoretical background
    print("This is a classic problem of stress concentration in mechanics of materials.")
    print("According to linear elasticity theory, the stress at a geometric discontinuity (like a notch or crack) is amplified.")
    print("For a sharp notch, the maximum stress at the tip (σ_max) is related to the nominal stress (σ_nominal) by a stress concentration factor, Kt.")
    
    # 3. Define the stress concentration factor and the maximum stress equation
    # The stress concentration factor Kt for a deep, sharp notch is often approximated as: Kt = 1 + 2 * sqrt(d / ρ)
    # where 'd' is the notch depth and 'ρ' is the notch tip radius.
    Kt = 1 + 2 * sp.sqrt(d / rho)
    sigma_max_expr = sigma_nom * Kt

    print(f"\nThe equation for the maximum stress at the notch tip is: σ_max = σ_nominal * (1 + 2 * sqrt(d / ρ))")

    # 4. Analyze the case of a perfectly sharp tip
    print("\nA 'perfectly sharp' wedge implies the radius of curvature at the tip 'ρ' approaches zero.")
    print("To find the theoretical stress, we calculate the limit of the σ_max expression as ρ -> 0.")

    # 5. Calculate the limit using sympy
    theoretical_stress = sp.limit(sigma_max_expr, rho, 0, dir='+')
    print(f"\nThe calculated limit of σ_max as ρ -> 0 is: {theoretical_stress}")
    
    print("\nConclusion: Theoretically, the stress at the tip of a perfectly sharp wedge (point A) is infinite.")

    # 6. Numerical Demonstration as requested by the user prompt
    print("\n--- Numerical Demonstration ---")
    # Assign some hypothetical values for demonstration
    hypothetical_sigma_nom = 100 # in MPa
    hypothetical_d = 5 # in mm
    hypothetical_rho = 0.0001 # a very small tip radius in mm

    print(f"Let's assume a nominal stress σ_nominal = {hypothetical_sigma_nom} MPa and a notch depth d = {hypothetical_d} mm.")
    print(f"For a very sharp tip radius ρ = {hypothetical_rho} mm, the stress concentration factor and maximum stress would be:")

    # Substitute numerical values
    kt_val = 1 + 2 * math.sqrt(hypothetical_d / hypothetical_rho)
    sigma_max_val = hypothetical_sigma_nom * kt_val
    
    # Print the final equation with each number, as requested.
    print("\nFinal Equation with example values:")
    print(f"Maximum Stress = Nominal Stress * (1 + 2 * sqrt(Depth / Radius))")
    print(f"{sigma_max_val:.2f} MPa = {hypothetical_sigma_nom} * (1 + 2 * sqrt({hypothetical_d} / {hypothetical_rho}))")
    print(f"{sigma_max_val:.2f} MPa = {hypothetical_sigma_nom} * ({kt_val:.2f})")
    print("\nAs the radius ρ gets even smaller, this value approaches infinity.")

if __name__ == '__main__':
    solve_stress_concentration()
