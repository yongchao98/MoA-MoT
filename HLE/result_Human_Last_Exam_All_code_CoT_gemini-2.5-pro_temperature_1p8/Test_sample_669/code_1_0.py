import math

def main():
    """
    Calculates the value of a_n,k,l mod p based on the given recurrence and parameters.
    """
    p = 21023

    # The recurrence relation leads to the generating function A_n(x,y) = P(x,y)^n, where
    # P(x,y) = 12 + 3x + 75y + 27x^2*y^2.
    # We need to find the coefficient of x^k * y^l in P(x,y)^n mod p.
    # By a generalization of Lucas's Theorem, this can be broken down based on the
    # base-p digits of n, k, and l.

    # The base-p digits (n_i, k_i, l_i) are periodic with period 3. We calculate
    # the corresponding coefficient C_i = [x^k_i y^l_i] P(x,y)^n_i for each pattern.

    # For i = 2, 5, 8,...: (n_i, k_i, l_i) = (2, 1, 1)
    # C_2 = [x^1 y^1] P(x,y)^2. This term comes from 2 * (3x) * (75y).
    c2_val = 2 * 3 * 75
    C2 = c2_val % p

    # For i = 1, 4, 7,...: (n_i, k_i, l_i) = (3, 1, 2)
    # C_1 = [x^1 y^2] P(x,y)^3. This term comes from (3!/(1!*2!)) * (3x)^1 * (75y)^2.
    c1_val = 3 * 3 * (75**2)
    C1 = c1_val % p

    # For i = 0, 3, 6,...: (n_i, k_i, l_i) = (5, 2, 2)
    # C_0 = [x^2 y^2] P(x,y)^5. This comes from two terms:
    # 1. From (5!/(1!2!2!)) * 12^1 * (3x)^2 * (75y)^2
    # 2. From (5!/(4!1!)) * 12^4 * (27x^2*y^2)^1
    c0_term1 = (math.factorial(5) // (math.factorial(1) * math.factorial(2) * math.factorial(2))) * 12 * (3**2) * (75**2)
    c0_term2 = (math.factorial(5) // (math.factorial(4) * math.factorial(1))) * (12**4) * 27
    c0_val = c0_term1 + c0_term2
    C0 = c0_val % p

    # The number of occurrences of each pattern up to the given limit means
    # a_n,k,l = C0^((p-1)/2+1) * C1^((p-1)/2+1) * C2^((p-1)/2) mod p
    # which simplifies to C0 * C1 * (Legendre symbol of (C0*C1*C2)) mod p.

    # Calculate the product for the Legendre symbol
    V = (C0 * C1 * C2) % p

    # Calculate Legendre symbol (V/p) = V^((p-1)/2) mod p
    # The result is 1 if V is a quadratic residue, p-1 (i.e., -1) if non-residue.
    L = pow(V, (p - 1) // 2, p)
    L_disp = 1 if L == 1 else -1

    # Final Result
    final_result = (C0 * C1 * L) % p
    
    # Print the equation as requested
    print("The coefficients for the digit patterns are:")
    print(f"C0 = [x^2 y^2]P(x,y)^5 = {c0_val} mod {p} = {C0}")
    print(f"C1 = [x^1 y^2]P(x,y)^3 = {c1_val} mod {p} = {C1}")
    print(f"C2 = [x^1 y^1]P(x,y)^2 = {c2_val} mod {p} = {C2}")
    print("\nThe final value is calculated as a_n,k,l = (C0 * C1 * (Legendre(C0*C1*C2, p))) mod p.")
    print("\nFinal equation:")
    print(f"({C0} * {C1} * {L_disp}) mod {p} = {final_result}")

if __name__ == "__main__":
    main()