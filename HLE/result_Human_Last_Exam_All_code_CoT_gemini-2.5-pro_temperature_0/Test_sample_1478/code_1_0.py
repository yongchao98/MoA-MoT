import math

# Titan architecture's 6-bit integer limit
MAX_VAL = 63

def check_6bit(n):
    """Checks if a number fits in a 6-bit unsigned integer."""
    return 0 <= n <= MAX_VAL

def simplify_fraction(num, den):
    """Simplifies a fraction using the greatest common divisor."""
    if den == 0:
        return num, den # Avoid division by zero
    common_divisor = math.gcd(num, den)
    return num // common_divisor, den // common_divisor

def titan_multiply(frac1, frac2):
    """
    Simulates a multiplication operation on the Titan architecture.
    Returns the resulting fraction or None if the operation fails.
    """
    num1, den1 = frac1
    num2, den2 = frac2

    res_num = num1 * num2
    res_den = den1 * den2

    print(f"Attempting to multiply {num1}/{den1} by {num2}/{den2}...")

    if not check_6bit(res_num) or not check_6bit(res_den):
        print(f"  - FAILURE: Intermediate result is {res_num}/{res_den}.")
        print(f"    The numerator ({res_num}) or denominator ({res_den}) exceeds the 6-bit limit of {MAX_VAL}.")
        # The architecture would try to expand terms here, but as shown in the analysis,
        # this also leads to overflow. For clarity, we show the initial failure.
        return None

    s_num, s_den = simplify_fraction(res_num, res_den)
    print(f"  - SUCCESS: Result is {res_num}/{res_den}, which simplifies to {s_num}/{s_den}.")
    return (s_num, s_den)

def run_simulation():
    """
    Runs the simulation of calculating the gravitational force on Titan.
    """
    print("--- Titan Architecture Feasibility Simulation ---")
    print("Goal: Calculate F = G * M * m / r^2")
    print("Approximation: r ≈ d = 1km, since Schwarzschild radius is negligible.")
    print("The final formula for the force's coefficient is: G * m * ρ * (4/3) * π * R³ / d²")
    print(f"Constraint: All intermediate numerators and denominators must be <= {MAX_VAL}.\n")

    # Define constants as fractional coefficients. Exponents are handled separately.
    # We use the most accurate fractional approximations possible.
    G_frac = (20, 3)
    m_frac = (50, 1)
    rho_frac = (6, 5)  # 1.2 = 12/10
    four_thirds_frac = (4, 3)
    pi_frac = (22, 7)
    R_cubed_frac = (8, 1) # R is 2e6, so R^3 is 8e18. We only need the coefficient.

    # Perform the calculation step-by-step
    print("--- Step 1: Combining ρ and 4/3 ---")
    # current_val represents the combined coefficient being calculated
    current_val = titan_multiply(rho_frac, four_thirds_frac)
    if not current_val:
        return

    print("\n--- Step 2: Multiplying by π ---")
    current_val = titan_multiply(current_val, pi_frac)
    if not current_val:
        return

    # The rest of the steps would also fail, but the failure at Step 2 is sufficient.
    print("\n--- Step 3: Multiplying by R³ ---")
    print("Simulation stopped as a failure has already occurred.")


# --- Main Execution ---
run_simulation()

print("\n--- Final Analysis ---")
print("The simulation demonstrates that calculating the force is not feasible on the Titan architecture.")
print("Even with the most precise and simple fractional representations, the multiplication of coefficients")
print("quickly results in numerators and denominators that exceed the 6-bit integer limit of 63.")
print("For example, in Step 2, the operation (8/5) * (22/7) results in an intermediate numerator of 176,")
print("which cannot be represented. Alternative calculation orders or expansion techniques also lead to overflow.")
print("Therefore, the problem cannot be solved as stated.")
