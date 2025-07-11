import math

def titan_multiply(val1, val2, exp):
    """
    Simulates Titan's multiplication, handling overflows by
    factoring out powers of 10 and approximating.
    Returns the new value and new exponent.
    """
    res = val1 * val2
    if res <= 63:
        print(f"* {val2} = {res} * 10^{exp}")
        return res, exp

    # Handle overflow
    original_res = res
    new_exp = exp
    while res > 63:
        # Find best factor of 10 to extract
        if res % 10 == 0:
             res //= 10
             new_exp += 1
        else: # Need to approximate
            # Approximate to the nearest 10
            approx_res = round(res / 10.0)
            res = approx_res
            new_exp += 1
    
    print(f"* {val2} = {original_res} -> approx {res} * 10^{new_exp}")
    return res, new_exp

def solve_pandora_gravity():
    """
    Calculates the gravitational force on a probe near Pandora-turned-black-hole
    using the Titan 6-bit architecture simulation.
    """
    # --- Step 1: Define constants and approximations ---
    # G = 6.674e-11 approx 20/3 e-11
    # m = 50 kg
    # r = 2000 km = 2e6 m => r^3 = 8e18 m^3
    # rho = 1200 kg/m^3 = 60 * 20 kg/m^3
    # d = 1 km = 1e3 m => d^2 = 1e6 m^2
    # pi approx 3/1
    
    # --- Step 2: Write the full equation ---
    print("The force F is calculated as F = G * m * M / d^2")
    print("where Mass M = Volume * Density = (4/3 * pi * r^3) * rho.")
    print("\nUsing 6-bit fractional approximations:")
    print("F = (G) * (m) * (4/3) * (pi) * (r^3) * (rho) / (d^2)")
    print("F = (20/3 * 10^-11) * (50) * (4/3) * (3) * (8 * 10^18) * (60 * 20) / (10^6)\n")

    # --- Step 3: Simplify the expression algebraically ---
    # Combine exponents: -11 + 18 - 6 = 1. So, * 10^1
    # Combine fractions: (20/3) * (50/1) * (4/3) * (3/1) * (60/1) * (20/1) * (8/1)
    # Denominators: 3 * 3 = 9
    # Numerators have factors 3 and 60. 3*60 = 180. 180 / 9 = 20.
    # All denominators cancel out.
    
    print("Simplifying the expression by cancelling terms:")
    print("The denominators (3*3=9) are cancelled by numerators (3*60=180).")
    print("The final exponent is 10^(-11+18-6) = 10^1.")
    print("The remaining numerator product is: 20 * 50 * 4 * 8 * 20 * 20\n")

    # --- Step 4: Simulate the sequential multiplication ---
    print("Calculating the numerator product with 6-bit constraints:")
    
    numerators = [20, 50, 4, 8, 20, 20]
    current_val = 1
    current_exp = 0
    print(f"Start with 1")

    # The first term is handled outside the loop to establish the start
    current_val = numerators.pop(0)
    print(f"-> {current_val} * 10^{current_exp}")
    
    for term in numerators:
        current_val, current_exp = titan_multiply(current_val, term, current_exp)

    final_num = current_val
    final_exp = current_exp + 1 # Add the exponent from the formula simplification
    final_den = 1 # Denominators all cancelled

    # --- Step 5: Final Result ---
    print("\nFinal Result:")
    print(f"The numerator calculation results in {final_num} * 10^{current_exp}.")
    print("Combining with the initial exponent of 10^1, the final force is:")
    print(f"F = {final_num}/{final_den} * 10^{final_exp} N")
    
    # --- Error calculation ---
    # True numerator value: 20*50*4*8*20*20 = 12,800,000 = 1.28e7
    # Approx numerator value: 13e6 = 1.3e7
    # Error from 128->130 approx: (1.3e7 - 1.28e7) / 1.28e7
    error_approx = (1.3 - 1.28) / 1.28
    
    # Error from pi approx: (math.pi - 3) / math.pi
    error_pi = (math.pi - 3) / math.pi
    
    total_error = error_approx + error_pi
    
    # The prompt expects an answer block at the end.
    # The format is Y[e] where e is the smallest relative error.
    answer = f"Y[{total_error*100:.1f}]"
    print(f"\n<<<__{answer}__>>>")


solve_pandora_gravity()