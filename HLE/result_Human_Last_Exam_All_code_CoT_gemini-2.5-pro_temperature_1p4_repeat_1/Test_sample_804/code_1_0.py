import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def titan_multiply(f1, f2):
    """
    Performs multiplication under Titan's 5-bit constraints.
    A fraction is represented as a tuple (numerator, denominator).
    Returns the resulting fraction or None if the operation is illegal.
    """
    n1, d1 = f1
    n2, d2 = f2

    # Check initial constraints
    if not all(0 <= x <= 31 for x in [n1, d1, n2, d2]):
        print(f"Error: Input fraction contains value > 31. Invalid inputs: {f1}, {f2}")
        return None
    
    # Perform cross-simplification as per the example
    g1 = gcd(n1, d2)
    n1, d2 = n1 // g1, d2 // g1
    
    g2 = gcd(n2, d1)
    n2, d1 = n2 // g2, d1 // g2
    
    # Calculate final numerator and denominator
    res_n = n1 * n2
    res_d = d1 * d2
    
    # Check if the result violates 5-bit constraints
    if res_n > 31 or res_d > 31:
        print(f"Failed to multiply {f1} by {f2}: Result ({res_n}, {res_d}) exceeds 5-bit limit.")
        return None
        
    return (res_n, res_d)

def demonstrate_calculation():
    """
    Attempts to calculate the gravitational force step-by-step
    to demonstrate the failure points.
    """
    print("Attempting to calculate F = (4/3) * pi * G * rho_shell * m_probe * r_total")
    print("All intermediate numerators and denominators must be <= 31.\n")
    
    # --- Define constants as Titan fractions ---
    four_thirds = (4, 3)
    m_probe = (30, 1)

    # These constants cannot be represented as 5-bit fractions due to their scale.
    # G = 6.674e-11, rho_shell = 300, r_total = 1,000,000
    # Their product is G * rho_shell * r_total ~= 0.02
    # The smallest non-zero representable fraction is 1/31 ~= 0.032.
    # So this term would be rounded to 0, making the force 0.
    print("Problem 1: Representing constants")
    print("G, rho_shell, and r_total cannot be represented due to their large/small scale.")
    print("Their product G*rho*r ~= 0.02, which is smaller than the smallest non-zero fraction 1/31.")
    print("This term would be rounded to 0, resulting in F=0.\n")

    # --- Attempt multiplication chain, showing failure ---
    print("Problem 2: Multiplying the representable constants")

    # Try with a good approximation for pi
    pi_accurate = (22, 7)
    print(f"Attempt 1: Using pi_T = {pi_accurate}")
    step1_res = titan_multiply(four_thirds, pi_accurate)
    if step1_res is None:
        print("Conclusion: Cannot multiply 4/3 by 22/7. This path fails.\n")

    # Try with a less accurate pi that allows multiplication
    pi_simplified = (3, 1)
    print(f"Attempt 2: Approximating pi_T = {pi_simplified} to enable the first step.")
    step1_res = titan_multiply(four_thirds, pi_simplified)
    if step1_res:
        print(f"Step 1 (4/3 * 3/1) successful. Result: {step1_res}")
        
        print("\nNow, multiplying by m_probe_T = (30, 1)")
        step2_res = titan_multiply(step1_res, m_probe)
        if step2_res is None:
            print("Conclusion: Cannot multiply the result by 30/1. This path also fails.")
            
    print("\nSince the calculation fails at fundamental steps due to the 5-bit constraints,")
    print("it is not possible to compute the gravitational force on Pandora using Titan.")

    # The problem asks for the equation, we show the one we tried to solve
    # F = 0 * 4/3 * 3/1 * 30/1, which is not right.
    # We state the calculation is impossible.
    print("\nFinal attempted equation:")
    print("F = (G * rho * r)_T * (4/3)_T * (pi)_T * (m_probe)_T")
    print("This equation cannot be resolved because:")
    print("1. The term (G*rho*r)_T would be 0.")
    print("2. The product of (4/3)_T * (pi)_T * (m_probe)_T exceeds the 5-bit limit.")
    print("\nFinal Answer: <<<N0>>>")


if __name__ == '__main__':
    demonstrate_calculation()