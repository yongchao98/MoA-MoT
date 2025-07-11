import math

def solve_magnetic_shell(R1, R2):
    """
    Calculates the required permeability (mu) of a cylindrical shell and the
    resulting interior magnetic field (H_int) for magnetic cloaking.

    This problem as stated only has a trivial solution (μ = μ_0). To find a
    non-trivial solution, it is assumed that the problem describes a
    magnetic cloak, which imposes a condition on the permeability of the
    interior region (μ_int), such that μ * μ_int = μ_0^2.

    The solution is valid only under the geometric condition R2 <= sqrt(2)*R1.

    Args:
        R1 (float): Inner radius of the cylindrical shell.
        R2 (float): Outer radius of the cylindrical shell.
    """
    
    # --- Check for valid inputs and physical constraints ---
    if not isinstance(R1, (int, float)) or not isinstance(R2, (int, float)) or R1 <= 0 or R2 <= 0:
        print("Error: Radii R1 and R2 must be positive numbers.")
        return
        
    if R1 >= R2:
        print(f"Error: Inner radius R1 ({R1}) must be smaller than outer radius R2 ({R2}).")
        return

    condition_val = 2 * R1**2 - R2**2
    if condition_val < 0:
        print(f"Error: No real solution for permeability exists for these radii (R1={R1}, R2={R2}).")
        print(f"The condition R2 <= sqrt(2)*R1 is not met. (sqrt(2)*R1 is approx {math.sqrt(2)*R1:.4f})")
        return

    # --- Calculations ---
    R1_sq = R1**2
    R2_sq = R2**2
    
    sqrt_term_val = math.sqrt(condition_val)
    numerator = R1_sq + R2 * sqrt_term_val
    denominator = R2_sq - R1_sq
    
    # k is the relative permeability μ/μ₀
    k = numerator / denominator

    # C is the factor relating the interior field H_int to the external field H_0
    C_numerator = k * (k + 1)
    C_denominator = k**2 + 1
    C = C_numerator / C_denominator

    # --- Output the results ---
    print(f"For a cylindrical shell with inner radius R1 = {R1} and outer radius R2 = {R2}:\n")

    print("1. Required Permeability of the Shell (μ):")
    print("------------------------------------------")
    print("The relative permeability k = μ/μ₀ is determined by the formula:")
    print("μ/μ₀ = [R₁² + R₂ * sqrt(2*R₁² - R₂²)] / (R₂² - R₁²)\n")
    print("Substituting the given values:")
    print(f"μ/μ₀ = [{R1_sq:.4f} + {R2:.4f} * sqrt(2*{R1_sq:.4f} - {R2_sq:.4f})] / ({R2_sq:.4f} - {R1_sq:.4f})")
    print(f"μ/μ₀ = [{R1_sq:.4f} + {R2:.4f} * sqrt({condition_val:.4f})] / {denominator:.4f}")
    print(f"μ/μ₀ = [{R1_sq:.4f} + {R2:.4f} * {sqrt_term_val:.4f}] / {denominator:.4f}")
    print(f"μ/μ₀ = [{R1_sq:.4f} + {R2 * sqrt_term_val:.4f}] / {denominator:.4f}")
    print(f"μ/μ₀ = {numerator:.4f} / {denominator:.4f}\n")
    
    print("Final result for permeability:")
    print(f">>> μ = {k:.4f} * μ₀\n")


    print("2. Magnetic Field in the Interior Region (H_int):")
    print("------------------------------------------------")
    print("The interior magnetic field H_int is uniform, given by H_int = C * H₀ * x_hat, where:")
    print("C = (μ/μ₀) * (μ/μ₀ + 1) / ((μ/μ₀)² + 1)\n")
    print("Substituting the calculated value for μ/μ₀:")
    print(f"C = {k:.4f} * ({k:.4f} + 1) / (({k:.4f})² + 1)")
    print(f"C = {k:.4f} * {k + 1:.4f} / ({k**2:.4f} + 1)")
    print(f"C = {C_numerator:.4f} / {C_denominator:.4f}\n")

    print("Final result for the interior field:")
    print(f">>> H_int = {C:.4f} * H₀ * x_hat")

# --- Define the radii for the problem ---
R1_val = 1.0
R2_val = 1.2

# --- Solve the problem and print the results ---
solve_magnetic_shell(R1_val, R2_val)
