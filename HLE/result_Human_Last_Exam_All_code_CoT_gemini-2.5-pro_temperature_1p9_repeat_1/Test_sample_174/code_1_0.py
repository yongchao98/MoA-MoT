from fractions import Fraction

def compute_euler_characteristic():
    """
    Computes the orbifold Euler characteristic of the moduli space of smooth plane quartics.
    """

    # Step 1: Compute e(M_3) using the Harer-Zagier formula.
    # The 6th Bernoulli number, B_6, is 1/42.
    B6 = Fraction(1, 42)
    
    # Zeta function at -5: zeta(-5) = -B_6 / 6
    zeta_minus_5 = -B6 / 6
    
    # Euler characteristic of M_3: e(M_3) = zeta(-5) / (-2)
    e_M3 = zeta_minus_5 / -2

    # Step 2: Use the known value for e(H_3).
    # From the literature, the Euler characteristic of the hyperelliptic locus H_3 is 1/45.
    e_H3 = Fraction(1, 45)

    # Step 3: Subtract to find the result.
    result = e_M3 - e_H3
    
    # --- Outputting the explanation and result ---
    
    print("The orbifold Euler characteristic of the quotient stack [U/G] is calculated as e(M_3) - e(H_3).")
    print("\nFirst, we calculate e(M_3):")
    print(f"  Bernoulli number B_6 = {B6}")
    print(f"  zeta(-5) = -B_6 / 6 = {zeta_minus_5}")
    print(f"  e(M_3) = zeta(-5) / -2 = {e_M3}")
    
    print("\nSecond, we use the known value for e(H_3):")
    print(f"  e(H_3) = {e_H3}")

    print("\nFinally, we compute the difference:")
    
    # Format the equation string.
    equation_str = f"e([U/G]) = {e_M3} - {e_H3} = {result}"
    
    # Replace Fraction() representation for a cleaner output.
    equation_str = equation_str.replace("Fraction(", "").replace(")", "")
    
    print("The final equation is:")
    print(f"  {e_M3} - {e_H3} = {result}")

compute_euler_characteristic()