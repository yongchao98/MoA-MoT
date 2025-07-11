from fractions import Fraction
import math

def rational_to_cf(p, q):
    """Computes the simple continued fraction of a rational number p/q."""
    coeffs = []
    if q == 0:
        return []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def evaluate_cf(coeffs):
    """Evaluates a finite simple continued fraction from a list of coefficients."""
    if not coeffs:
        return Fraction(0)
    val = Fraction(coeffs[-1])
    for i in range(len(coeffs) - 2, -1, -1):
        val = coeffs[i] + 1 / val
    return val

def cf_quadratic_irrational(D, P0, Q0):
    """
    Computes the continued fraction of a quadratic irrational of the form (P0 + sqrt(D)) / Q0.
    The algorithm requires that Q0 divides (D - P0^2). This function handles the case
    where this is not initially true by transforming the number.
    """
    if (D - P0**2) % Q0 != 0:
        P0_new = Q0 * P0
        D_new = D * Q0**2
        Q0_new = Q0 * Q0
        P0, D, Q0 = P0_new, D_new, Q0_new

    coeffs = []
    visited = {}
    Pk, Qk = P0, Q0
    
    while True:
        # If we see the same (P, Q) pair, we have found the repeating cycle.
        if (Pk, Qk) in visited:
            start_index = visited[(Pk, Qk)]
            non_repeating = coeffs[:start_index]
            repeating = coeffs[start_index:]
            return non_repeating, repeating

        visited[(Pk, Qk)] = len(coeffs)
        
        # Calculate the next coefficient
        ak = int((Pk + math.sqrt(D)) // Qk)
        coeffs.append(ak)
        
        # Calculate the next P and Q values
        P_next = ak * Qk - Pk
        Q_next = (D - P_next**2) // Qk
        
        Pk, Qk = P_next, Q_next

def main():
    """Main function to compute and print the result."""
    # Step 1: Define the rational number
    p, q = 4, 7
    print(f"Computing the continued fraction for the generalized Markov number m_{{{p}/{q}}}.")

    # Step 2: Get the continued fraction of p/q = [0; a1, a2, ...]
    cf_coeffs = rational_to_cf(p, q)[1:] # We need the part after the initial zero
    n = len(cf_coeffs)

    # Step 3: Calculate the Markov Constant C_{p/q}
    max_C = Fraction(0)
    for k in range(1, n + 1):
        lambda_coeffs = cf_coeffs[k-1:]
        lambda_val = evaluate_cf(lambda_coeffs)
        
        mu_coeffs = cf_coeffs[:k-1]
        mu_coeffs.reverse()
        if not mu_coeffs:
            mu_val = Fraction(0)
        else:
            mu_val = 1 / evaluate_cf(mu_coeffs)
        
        current_C = lambda_val + mu_val
        if current_C > max_C:
            max_C = current_C

    # Step 4: Calculate m_{p/q} = sqrt(C^2 + 4)
    C_num = max_C.numerator
    C_den = max_C.denominator
    D_m = C_num**2 + 4 * C_den**2
    P_m = 0
    Q_m = C_den
    
    # Step 5: Compute the continued fraction of m_{p/q}
    non_rep, rep = cf_quadratic_irrational(D_m, P_m, Q_m)

    # Step 6: Print the result clearly
    print(f"\nThe generalized Markov number m_{p}/{q} is sqrt({D_m})/{Q_m}.")
    print("Its continued fraction has the form [a_0; (a_1, a_2, ...)], where the parentheses indicate the repeating part.")
    
    # Output each number in the final equation
    print("\nThe numbers in the final continued fraction are:")
    print(f"The non-repeating part is a_0 = {non_rep[0]}.")
    rep_str = ", ".join(map(str, rep))
    print(f"The repeating part is ({rep_str}).")
    
    full_cf_str = f"[{non_rep[0]}; ({rep_str})]"
    print(f"\nThus, the full continued fraction is: {full_cf_str}")

if __name__ == '__main__':
    main()
