import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def simplify_fraction(num, den):
    """Simplifies a fraction by dividing both parts by their GCD."""
    if den == 0:
        raise ValueError("Denominator cannot be zero.")
    common = gcd(num, den)
    return num // common, den // common

def titan_mul(frac1, frac2, current_exponent):
    """
    Simulates Titan's 6-bit multiplication with reduction.
    frac1 and frac2 are tuples (numerator, denominator).
    Returns the new fraction and the updated exponent.
    """
    num1, den1 = frac1
    num2, den2 = frac2
    
    # Standard fraction multiplication
    new_num = num1 * num2
    new_den = den1 * den2
    
    # Titan's reduction rule for numbers > 63
    while new_num > 63:
        new_num //= 10
        current_exponent += 1
        
    while new_den > 63:
        new_den //= 10
        # Dividing the denominator by 10 is equivalent to multiplying the number by 10
        current_exponent += 1
        
    # Simplify the resulting fraction
    return simplify_fraction(new_num, new_den), current_exponent

def main():
    """
    Main function to perform the Titan calculation for the gravitational force.
    """
    # Step 1: Define all constants and variables in Titan's format (fraction, exponent)
    # F = G * (rho * 4/3 * pi * R^3) * m_probe / r^2
    
    # Gravitational Constant G ≈ 6.67e-11 N m^2/kg^2 -> use 20/3 * 10^-11
    G = ((20, 3), -11)
    
    # Density rho = 1200 kg/m^3 -> use 12/1 * 10^2
    rho = ((12, 1), 2)
    
    # Constant 4/3
    four_thirds = ((4, 3), 0)
    
    # Pi ≈ 22/7
    pi = ((22, 7), 0)
    
    # Radius R = 2000 km = 2e6 m. We need R^3 = (2e6)^3 = 8e18 m^3
    R_cubed = ((8, 1), 18)
    
    # Probe mass m = 50 kg
    m_probe = ((50, 1), 0)
    
    # Distance r ≈ 1 km = 1000 m. We need r^2 = (1000)^2 = 1e6 m^2
    r_squared = ((1, 1), 6)
    
    # List of all terms to be multiplied in the numerator
    terms = [G, rho, four_thirds, pi, R_cubed, m_probe]
    
    # Step 2: Calculate the initial exponent
    # Sum of exponents from terms, minus exponent from denominator (r_squared)
    final_exponent = sum(term[1] for term in terms) - r_squared[1]
    
    # Step 3: Sequentially multiply the mantissas, applying Titan's rules
    print("Starting Titan Feasibility Calculation...\n")
    print("Objective: Calculate Force F = G * rho * (4/3) * pi * R^3 * m_probe / r^2")
    
    # Initialize the running product with the first term
    running_mantissa = terms[0][0]
    print(f"MOV AX, {running_mantissa[0]}/{running_mantissa[1]}")
    
    # Loop through the rest of the terms
    for i in range(1, len(terms)):
        term_mantissa = terms[i][0]
        print(f"MUL AX, {term_mantissa[0]}/{term_mantissa[1]}")
        running_mantissa, final_exponent = titan_mul(running_mantissa, term_mantissa, final_exponent)
        print(f"  Intermediate Result: AX = {running_mantissa[0]}/{running_mantissa[1]}, Exponent = {final_exponent}")
    
    # The final operation is division by r_squared, which is 1/1, so it doesn't change the result.
    # We can represent division by A as multiplication by 1/A.
    r_squared_inv = (r_squared[0][1], r_squared[0][0]) # (1,1)
    print(f"DIV AX, {r_squared[0][0]}/{r_squared[0][1]} (equivalent to MUL AX, {r_squared_inv[0]}/{r_squared_inv[1]})")
    running_mantissa, final_exponent = titan_mul(running_mantissa, r_squared_inv, final_exponent)
    print(f"  Intermediate Result: AX = {running_mantissa[0]}/{running_mantissa[1]}, Exponent = {final_exponent}")
    
    final_num, final_den = running_mantissa
    
    print("\n--- RED AX ---")
    print("Final Titan Calculation Result:")
    print(f"F = ({G[0][0]}/{G[0][1]}) * ({rho[0][0]}/{rho[0][1]}) * ({four_thirds[0][0]}/{four_thirds[0][1]}) * ({pi[0][0]}/{pi[0][1]}) * ({R_cubed[0][0]}/{R_cubed[0][1]}) * ({m_probe[0][0]}/{m_probe[0][1]}) / ({r_squared[0][0]}/{r_squared[0][1]}) * 10^({G[1]} + {rho[1]} + {four_thirds[1]} + {pi[1]} + {R_cubed[1]} + {m_probe[1]} - {r_squared[1]})")
    print(f"\nSimplified Force (F) = {final_num}/{final_den} x 10^{final_exponent} N\n")

    # Step 4: Calculate error
    # "True" value calculation with higher precision
    G_true = 6.674e-11
    rho_true = 1.2e3
    R_true = 2e6
    m_probe_true = 50
    r_true = 1e3
    M_true = rho_true * (4/3) * math.pi * (R_true**3)
    F_true = (G_true * M_true * m_probe_true) / (r_true**2)

    # Calculated value from Titan simulation
    F_calc = (final_num / final_den) * (10**final_exponent)
    
    # Relative error
    relative_error = abs(F_calc - F_true) / F_true
    error_percent = round(relative_error * 100, 1)

    print("--- Verification ---")
    print(f"Titan Calculated Value: {F_calc:.4e} N")
    print(f"High-Precision Value: {F_true:.4e} N")
    print(f"Smallest Relative Error Achieved: {error_percent}%")

if __name__ == '__main__':
    main()
