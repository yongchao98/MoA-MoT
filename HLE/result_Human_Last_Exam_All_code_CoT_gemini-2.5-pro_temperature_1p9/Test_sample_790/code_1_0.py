import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def multiply_fractions(f1, f2):
    """
    Multiplies two fractions, simplifying before multiplying to adhere to Titan constraints.
    Returns the resulting fraction if valid, otherwise None.
    """
    n1, d1 = f1
    n2, d2 = f2

    # Simplify before multiplying to avoid overflow
    common_divisor1 = gcd(n1, d2)
    n1_s = n1 // common_divisor1
    d2_s = d2 // common_divisor1

    common_divisor2 = gcd(n2, d1)
    n2_s = n2 // common_divisor2
    d1_s = d1 // common_divisor2

    num = n1_s * n2_s
    den = d1_s * d2_s
    
    # Check if the result violates the 5-bit constraint
    if num > 31 or den > 31:
        return None
        
    return (num, den)

def format_fraction(f):
    """Formats a fraction tuple as a string 'n/d'."""
    return f"{f[0]}/{f[1]}"

def solve_titan_problem():
    """
    Solves the problem using Titan computer architecture constraints.
    """
    # 1. Define true values for calculating error
    r_true = 0.005  # radius in meters
    rho_true = 900000 # density in kg/m^3 (0.9 kg/cm^3)
    pi_true = math.pi
    m_true = (4/3) * pi_true * (r_true**3) * rho_true
    g_true = 9.8
    sqrt2_true = math.sqrt(2)
    F_true = 2 * sqrt2_true * m_true * g_true
    
    # 2. Choose workable Titan fraction approximations for the constants
    # m_true is ~0.471 kg. We approximate it as 1/2.
    m_f = (1, 2)
    # g_true is 9.8. We approximate it as 10/1.
    g_f = (10, 1)
    # sqrt2_true is ~1.414. We use the 7/5 approximation for calculability.
    sqrt2_f = (7, 5)
    # The initial factor of 2 is represented as 2/1.
    two_f = (2, 1)

    # 3. Perform the calculation with an order that keeps intermediate results valid.
    # Equation: F = 2 * sqrt(2) * m * g
    # We will calculate (2 * m) * (g * sqrt(2))
    
    # Calculate term1 = 2 * m
    term1 = multiply_fractions(two_f, m_f)
    if term1 is None:
        print("Calculation failed at term1. This combination is not viable.")
        return

    # Calculate term2 = g * sqrt(2)
    term2 = multiply_fractions(g_f, sqrt2_f)
    if term2 is None:
        print("Calculation failed at term2. This combination is not viable.")
        return

    # Calculate final result F = term1 * term2
    F_f = multiply_fractions(term1, term2)
    if F_f is None:
        print("Final calculation failed. This combination is not viable.")
        return
        
    F_titan_val = F_f[0] / F_f[1]
    
    # 4. Print the final equation as required
    print("Titan Calculation:")
    # The order shown here is for readability; the calculation was performed as above.
    expression = f"{format_fraction(two_f)} * {format_fraction(sqrt2_f)} * {format_fraction(m_f)} * {format_fraction(g_f)}"
    result_str = f"{format_fraction(F_f)}"
    print(f"{expression} = {result_str}")
    
    # 5. Calculate and print the error
    error = abs(F_titan_val - F_true)
    
    print("\nSummary:")
    print(f"Approximated Force (Titan) = {F_titan_val:.3f} N")
    print(f"More Precise Force = {F_true:.3f} N")
    print(f"Absolute Error = {error:.3f} N")
    
    # 6. Format final answer
    final_answer = f"Y[{error:.3f}]"
    print(f"\nFinal Answer: {final_answer}")
    print(f"\n<<<{final_answer}>>>")

solve_titan_problem()