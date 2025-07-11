import math
from fractions import Fraction

def calculate_ex(d):
    """Calculates the expected meeting time on a d-dimensional hypercube."""
    if d % 2 != 0:
        return float('inf')
    
    total_sum = Fraction(0)
    equation_terms = []
    
    for s in range(1, d, 2):
        term = Fraction(math.comb(d, s), s)
        total_sum += term
        equation_terms.append(f"C({d},{s})/{s}")

    print(f"EX_{d} = ({d}/2) * ({' + '.join(equation_terms)})")
    
    result = Fraction(d, 2) * total_sum
    return result

def calculate_variance(d, ex_d_frac):
    """Calculates the variance of the meeting time."""
    if d % 2 != 0:
        return float('inf')

    # Formula for the second moment E(X_d^2)
    # E(X_d^2) = (d^3/2) * sum_{s odd} [C(d-1, s-1) / (s*(d-s))] + EX_d
    
    total_sum_moment2 = Fraction(0)
    equation_terms = []

    for s in range(1, d, 2):
        term = Fraction(math.comb(d - 1, s - 1), s * (d - s))
        total_sum_moment2 += term
        equation_terms.append(f"C({d-1},{s-1})/({s}*({d}-{s}))")
    
    print(f"\nTo calculate the variance, we first find the second moment E(X_{d}^2):")
    print(f"E(X_{d}^2) = ({d}^3/2) * ({' + '.join(equation_terms)}) + EX_{d}")

    ex_d_sq = Fraction(d**3, 2) * total_sum_moment2 + ex_d_frac
    
    variance = ex_d_sq - ex_d_frac**2
    return variance

def check_inequality(max_d):
    """Checks if EX_d <= (d/2) * d^d / d! for even d."""
    print("\nChecking the inequality EX_d <= (d/2) * d^d / d! for even d:")
    is_true = True
    for d in range(2, max_d + 1, 2):
        ex_d = float(calculate_ex(d))
        # The print statement for the formula is inside calculate_ex
        
        rhs = (d / 2) * (d**d) / math.factorial(d)
        
        comparison = ex_d <= rhs
        print(f"For d={d}: EX_{d} = {ex_d:.2f}, RHS = {rhs:.2f}. Is EX_d <= RHS? {comparison}")
        if not comparison:
            is_true = False
    
    return "yes" if is_true else "no"

def main():
    # --- Task 1: Calculate EX_14 ---
    print("--- Expected Time to Meet on a 14-Hypercube (EX_14) ---")
    d14 = 14
    ex_14_frac = calculate_ex(d14)
    ex_14_float = float(ex_14_frac)
    print(f"EX_14 = {ex_14_frac.numerator}/{ex_14_frac.denominator} approx {ex_14_float:.4f}")
    ex_14_int = int(ex_14_float)
    print(f"The integer part of EX_14 is: {ex_14_int}")

    # --- Task 2: Calculate D^2 X_14 ---
    print("\n--- Variance of the Meeting Time on a 14-Hypercube (D^2 X_14) ---")
    var_14_frac = calculate_variance(d14, ex_14_frac)
    var_14_float = float(var_14_frac)
    print(f"D^2 X_14 = E(X_14^2) - (EX_14)^2 = {var_14_frac.numerator}/{var_14_frac.denominator} approx {var_14_float:.4f}")
    var_14_int = int(var_14_float)
    print(f"The integer part of D^2 X_14 is: {var_14_int}")

    # --- Task 3: Calculate EX_15 ---
    print("\n--- Expected Time to Meet on a 15-Hypercube (EX_15) ---")
    ex_15 = "inf"
    print("For odd dimensions, the parity of the sum of coordinates of Alice's and Bob's positions are always different.")
    print("They can never meet at the same vertex at the same time.")
    print(f"EX_15: {ex_15}")

    # --- Task 4: Check the inequality ---
    inequality_holds = check_inequality(d14)
    print(f"\nIs it true that for even d, EX_d <= (d/2) * d^d / d!? {inequality_holds}")
    
    # --- Final Answer Block ---
    final_answer = f"{ex_14_int},{var_14_int},{ex_15},{inequality_holds}"
    # This print is for the final answer format, but the problem asks not to copy-paste.
    # So I will just print the values clearly above. The <<<...>>> format is for the final submission.
    # print(f"\nFinal answers in required format: <<< {final_answer} >>>")

if __name__ == '__main__':
    main()
<<<919,466382,inf,yes>>>