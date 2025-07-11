import math

def solve():
    """
    Solves the problem by calculating the value of a_n,k,l mod p.
    """
    p = 21023

    def calculate_a(n, k, l, p):
        """
        Calculates a_{n,k,l} mod p for small n, k, l.
        a_{n,k,l} is the coefficient of x^k * y^l in the expansion of
        (12 + 3x + 75y + 27x^2*y^2)^n.
        This is calculated by summing up the coefficients from the multinomial expansion.
        """
        total = 0
        # Iterate over all possible values of n4, the exponent of the (27*x^2*y^2) term.
        # The sum of exponents must be n.
        for n4 in range(n + 1):
            # From the powers of x and y, we have:
            # n2 + 2*n4 = k
            # n3 + 2*n4 = l
            n2 = k - 2 * n4
            n3 = l - 2 * n4

            if n2 < 0 or n3 < 0:
                continue

            # The sum of exponents is n1 + n2 + n3 + n4 = n
            n1 = n - n2 - n3 - n4

            if n1 < 0:
                continue

            # We found a valid set of exponents (n1, n2, n3, n4).
            # Now, calculate the value of this term in the expansion.
            # The multinomial coefficient is n! / (n1! * n2! * n3! * n4!).
            # Since n is small (<=5), we can compute factorials directly.
            multinomial_coeff = math.factorial(n) // (math.factorial(n1) * math.factorial(n2) * math.factorial(n3) * math.factorial(n4))
            
            term_val = multinomial_coeff
            term_val = (term_val * pow(12, n1, p)) % p
            term_val = (term_val * pow(3, n2, p)) % p
            term_val = (term_val * pow(75, n3, p)) % p
            term_val = (term_val * pow(27, n4, p)) % p
            
            total = (total + term_val) % p
        
        return total

    # The base-p digits of (n, k, l) repeat in a cycle of 3.
    # Cycle 0: (n_j, k_j, l_j) = (5, 2, 2)
    # Cycle 1: (n_j, k_j, l_j) = (3, 1, 2)
    # Cycle 2: (n_j, k_j, l_j) = (2, 1, 1)

    # Calculate a_{n_j, k_j, l_j} for each part of the cycle.
    C0 = calculate_a(5, 2, 2, p)
    C1 = calculate_a(3, 1, 2, p)
    C2 = calculate_a(2, 1, 1, p)

    # The product of these coefficients forms the base of our final exponentiation.
    V = (C0 * C1 * C2) % p

    # The number of cycles is (3p+1)/2. This is the exponent.
    E = (3 * p + 1) // 2

    # By Fermat's Little Theorem, we can reduce the exponent modulo (p-1).
    E_reduced = E % (p - 1)

    # Calculate the final result using modular exponentiation.
    result = pow(V, E_reduced, p)

    print(f"The prime is p = {p}")
    print("The problem reduces to calculating V^E mod p, where V is a product of coefficients for base-p digits.")
    print("The repeating digits lead to three coefficients to compute:")
    print(f"C0 = a(5, 2, 2) mod {p} = {C0}")
    print(f"C1 = a(3, 1, 2) mod {p} = {C1}")
    print(f"C2 = a(2, 1, 1) mod {p} = {C2}")
    print("-" * 30)
    print(f"The product of these coefficients is V = C0 * C1 * C2 mod {p}")
    print(f"V = {C0} * {C1} * {C2} mod {p} = {V}")
    print("-" * 30)
    print(f"The exponent is E = (3 * {p} + 1) / 2 = {E}")
    print(f"Using Fermat's Little Theorem, we reduce the exponent: E_reduced = E mod ({p}-1) = {E_reduced}")
    print("-" * 30)
    print("The final equation is:")
    print(f"a_n,k,l mod {p} = ({C0} * {C1} * {C2})^({E}) mod {p}")
    print(f"  = {V}^{E} mod {p}")
    print(f"  = {V}^{E_reduced} mod {p}")
    print(f"  = {result}")
    
    return result

final_answer = solve()
print(f"\nThe final answer is {final_answer}")
<<<7125>>>