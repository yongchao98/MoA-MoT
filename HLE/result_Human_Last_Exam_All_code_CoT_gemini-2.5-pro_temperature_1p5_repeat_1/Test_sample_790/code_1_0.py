import math

# Helper function to represent fractions and simplify them
def simplify_fraction(n, d):
    """Simplifies a fraction n/d."""
    if d == 0:
        raise ValueError("Denominator cannot be zero")
    common_divisor = math.gcd(n, d)
    return (n // common_divisor, d // common_divisor)

def multiply_fractions(f1, f2, limit=31):
    """Multiplies two fractions (n1, d1) and (n2, d2) and checks constraints."""
    n1, d1 = f1
    n2, d2 = f2
    
    # Pre-simplification (as done in Titan)
    n1, d2_s = simplify_fraction(n1, d2)
    n2, d1_s = simplify_fraction(n2, d1)

    new_n = n1 * n2
    new_d = d1_s * d2_s
    
    if new_n > limit or new_d > limit:
        # This function signals that a manual approximation step is needed
        return None 
    
    return (new_n, new_d)

def print_step(message, fraction=None, value=None):
    if fraction:
        n, d = fraction
        print(f"{message} {(n, d)}")
    elif value is not None:
         print(f"{message} {value}")
    else:
        print(message)

# --- Main Calculation ---
def solve_titan_problem():
    # Target value for Force F, calculated with high precision for error calculation
    real_m = 0.9 * (4/3) * math.pi * (0.5**3) # mass in kg
    real_g = 9.8
    real_sqrt2 = math.sqrt(2)
    real_F = 2 * real_sqrt2 * real_m * real_g # ~13.08 N

    print("Step 1: Define initial constants as fractions (n, d) with n, d <= 31")
    # Constants for mass calculation
    rho = (9, 10)  # 0.9 kg/cm^3
    four_thirds = (4, 3)
    # Pi must be chosen carefully. A numerator > 10 makes m calculation fail.
    # π ≈ 3 is represented as (3, 1)
    pi = (3, 1) 
    r_cubed = (1, 8) # r = 0.5cm, r^3 = 0.125 = 1/8

    # Constants for force calculation
    g = (10, 1) # g ≈ 10 m/s^2
    two = (2, 1)
    # We will start with a common approximation for sqrt(2) ≈ 1.4
    sqrt2_approx = (7, 5)

    print_step("Density ρ:", rho)
    print_step("Approximation for 4/3:", four_thirds)
    print_step("Approximation for π:", pi)
    print_step("Radius cubed (r³):", r_cubed)
    print_step("Approximation for g:", g)
    print_step("Approximation for √2:", sqrt2_approx)
    print("-" * 30)

    print("Step 2: Calculate mass 'm' using fractional arithmetic.")
    # m = ρ * (4/3) * π * r³
    # (9/10) * (4/3) -> simplify 9/3 and 4/10 -> (3/5) * (2/1) = (6, 5)
    temp_m1 = multiply_fractions(rho, four_thirds)
    print_step("Intermediate m = ρ * (4/3) =", temp_m1)
    # (6/5) * (3/1) = (18, 5)
    temp_m2 = multiply_fractions(temp_m1, pi)
    print_step("Intermediate m = m * π =", temp_m2)
    # (18/5) * (1/8) -> simplify 18/8 -> (9/5) * (1/4) = (9, 20)
    m = multiply_fractions(temp_m2, r_cubed)
    print_step("Final calculated mass m =", m)
    print("-" * 30)

    print("Step 3: Calculate force F = 2 * m * g * √2")
    # To keep numerators small, we reorder: F = (2 * m) * g * sqrt(2)
    # Let's try an order that simplifies well: F = ( (2 * m) * g ) * √2
    temp_f1 = multiply_fractions(two, m)
    print_step("Intermediate F = 2 * m =", temp_f1, "(simplified from (18, 20))")
    
    temp_f2 = multiply_fractions(temp_f1, g)
    print_step("Intermediate F = (2 * m) * g =", temp_f2, "(simplified from (90, 10))")
    
    final_f_n, final_f_d = temp_f2
    sqrt_n, sqrt_d = sqrt2_approx
    
    print(f"\nStep 4: Constraint Maintenance")
    print(f"The next operation is ({final_f_n}, {final_f_d}) * ({sqrt_n}, {sqrt_d}).")
    print(f"The resulting numerator would be {final_f_n} * {sqrt_n} = {final_f_n * sqrt_n}, which is greater than 31.")
    print("This violates the 5-bit constraint.")
    
    print("\nTo proceed, we must replace an operand with a less precise fraction.")
    print(f"We need to approximate √2 ≈ 1.414 with a fraction (a/b) where {final_f_n}*a <= 31, so a <= {31//final_f_n}.")
    new_sqrt2_approx = (3, 2) # 1.5 is a good approximation and 9*3=27<=31
    print(f"We replace the √2 approximation (7, 5) with (3, 2), which represents 1.5.")
    
    final_F = multiply_fractions(temp_f2, new_sqrt2_approx)
    
    print("\nStep 5: Final Calculation")
    print_step(f"Final F = ({final_f_n}/{final_f_d}) * ({new_sqrt2_approx[0]}/{new_sqrt2_approx[1]}) =", final_F)
    
    final_val = final_F[0] / final_F[1]
    error = abs(final_val - real_F)
    
    print("-" * 30)
    print("Final Result:")
    print(f"The calculated force is ({final_F[0]}) / ({final_F[1]}) = {final_val:.1f} N")
    
    print(f"The absolute error is |{final_val:.3f} - {real_F:.3f}| = {error:.3f}")
    
    # Final answer in specified format
    final_answer = f"Y[{error:.3f}]"
    print(f"\nFinal Answer: {final_answer}")
    print(f"\n<<<Y[{error:.3f}]>>>")

solve_titan_problem()