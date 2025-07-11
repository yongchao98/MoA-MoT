import math

def calculate_norm_T_for_odd_n(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with an odd integer n.

    The function demonstrates the calculation via direct summation and verifies it
    with the simplified analytical formula 2**(n+1) - 1.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: Input 'n' must be a positive odd integer.")
        return

    print(f"Calculating the 1-norm of the correlation matrix T for n = {n}.")
    print("The formula is: ||T||_1 = (1 / (1 + 3^n)) * Σ_{w=1 to n+1} [ C(n+1, w) * |3^w + 3^n * (-1)^w| ]\n")

    # Numerator calculation via summation
    numerator_sum = 0
    n_plus_1 = n + 1
    
    print("Calculating the sum in the numerator term by term:")
    for w in range(1, n_plus_1 + 1):
        # Binomial coefficient C(n+1, w)
        comb = math.comb(n_plus_1, w)
        
        # Absolute value term |3^w + 3^n * (-1)^w|
        term_in_abs = 3**w + (3**n) * ((-1)**w)
        abs_val = abs(term_in_abs)
        
        term = comb * abs_val
        numerator_sum += term
        
        sign = "-" if w % 2 != 0 else "+"
        print(f"w={w:2d}: C({n_plus_1}, {w}) * |3^{w} {sign} 3^{n}| = {comb:<4} * {abs_val:<10} = {term}")

    # Denominator
    denominator = 1 + 3**n
    
    print(f"\nNumerator (Sum) = {numerator_sum}")
    print(f"Denominator (1 + 3^n) = 1 + {3**n} = {denominator}")
    
    # Final result from the sum
    final_result_sum = numerator_sum / denominator
    print(f"\nResult from sum: ||T||_1 = {numerator_sum} / {denominator} = {final_result_sum:.1f}\n")
    
    print("-" * 40)
    
    # Verification with the simplified formula
    print("Verifying with the simplified formula: 2^(n+1) - 1")
    final_result_formula = 2**(n_plus_1) - 1
    print(f"2^({n}+1) - 1 = 2^{n_plus_1} - 1 = {2**n_plus_1} - 1 = {final_result_formula}")

    if math.isclose(final_result_sum, final_result_formula):
        print("\nThe results from both methods match.")
    else:
        print("\nWarning: The results do not match. Please check the derivation.")

if __name__ == '__main__':
    # You can change this value to any positive odd integer
    odd_n = 5
    calculate_norm_T_for_odd_n(odd_n)