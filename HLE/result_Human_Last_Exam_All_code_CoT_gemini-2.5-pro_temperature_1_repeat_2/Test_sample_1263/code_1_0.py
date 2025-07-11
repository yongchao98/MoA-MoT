import math

def check_overflow(n, d):
    """Checks if a fraction n/d overflows the 4-bit constraint."""
    if n > 15 or d > 15:
        return True
    return False

def titan_multiply(f1_n, f1_d, f2_n, f2_d):
    """
    Simulates a multiplication on Titan.
    f1 and f2 are the two fractions (n/d).
    """
    res_n = f1_n * f2_n
    res_d = f1_d * f2_d
    
    # Titan rule: "immediately simplified"
    common_divisor = math.gcd(res_n, res_d)
    final_n = res_n // common_divisor
    final_d = res_d // common_divisor
    
    print(f"Multiplying {f1_n}/{f1_d} by {f2_n}/{f2_d}...")
    print(f"Intermediate result: {res_n}/{res_d}")
    print(f"Simplified result: {final_n}/{final_d}")
    
    if check_overflow(final_n, final_d):
        print("Result: OVERFLOW. Calculation cannot proceed.")
        return False
    else:
        print("Result: OK.")
        return True

def main():
    print("Analyzing feasibility of calculating Pandora's escape velocity on Titan.")
    print("-" * 20)
    
    # We need to calculate v_e^2, which involves the product of several constants.
    # The core numerical part is ~ 6.67 * 3.14 * 3.2
    
    # Let's use the best possible 4-bit fractional approximations for the mantissas.
    # G's mantissa (6.67) is best approximated by 13/2 = 6.5
    G_mantissa_n, G_mantissa_d = 13, 2
    
    # pi's mantissa (3.14) is best approximated by 13/4 = 3.25
    pi_mantissa_n, pi_mantissa_d = 13, 4
    
    print("Step 1: Multiply approximation of G by approximation of pi.")
    # The first step of the calculation would be to multiply G by pi.
    # Let's see if this is possible.
    if not titan_multiply(G_mantissa_n, G_mantissa_d, pi_mantissa_n, pi_mantissa_d):
        print("\nConclusion: The very first multiplication of the core constants causes an overflow.")
        print("The Titan architecture cannot handle numbers of this magnitude.")
        print("Therefore, it is not possible to calculate the escape velocity.")
        print("\nFinal Answer: N0")

main()