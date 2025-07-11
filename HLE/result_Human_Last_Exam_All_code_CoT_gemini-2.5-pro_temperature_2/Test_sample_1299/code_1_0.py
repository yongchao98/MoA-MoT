import math

def get_prime_factorization(n):
    """
    Calculates the prime factorization of a number n.
    Returns a dictionary mapping each prime factor to its exponent.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_and_explain(ell):
    """
    Calculates the cardinality of T_ell and U_ell, and prints the formulas.
    """
    print(f"--- Calculating for l = {ell} ---")

    # Part A and B are expressions for the same value, as |U_ell| = |T_ell|.
    # Part A asks for an expression in terms of l.
    # The derived formula is: |U_ell| = 1 if l=1, and d(l^2) - 1 if l > 1.
    # d(n) is the divisor function. This can be written as d(l^2) - 1 + delta_{l,1},
    # where delta is the Kronecker delta.

    # Part B asks for an expression in terms of the prime factorization l = p_1^e_1 * ... * p_s^e_s.
    # The derived formula is: |T_ell| = ( (2*e_1+1)...(2*e_s+1) ) - 1 + delta_{s,0},
    # where s is the number of distinct prime factors and delta_{s,0} is 1 if s=0 (l=1) and 0 otherwise.

    # Now we demonstrate the calculation for the given l.
    if ell == 1:
        s = 0
        exponents = []
        result = 1
        print("For l = 1, the number of distinct prime factors s = 0.")
        print("The formula is (empty product) - 1 + 1 = 1 - 1 + 1 = 1.")
    else:
        prime_factors = get_prime_factorization(ell)
        exponents = list(prime_factors.values())
        s = len(exponents)

        print(f"Prime factorization of {ell} is: {' * '.join([f'{p}^{e}' for p, e in prime_factors.items()])}")
        print(f"The exponents are: {exponents}")
        
        # Calculate d(ell^2) from exponents
        d_ell_sq_terms = [2 * e + 1 for e in exponents]
        
        # Build the equation string
        calc_str = f"|T_{ell}| = ({' * '.join(map(str, d_ell_sq_terms))}) - 1"
        
        d_ell_sq = 1
        for term in d_ell_sq_terms:
            d_ell_sq *= term
        
        result = d_ell_sq - 1
        
        if len(d_ell_sq_terms) > 1:
            calc_str += f" = {d_ell_sq} - 1"

        calc_str += f" = {result}"
        print("Calculation using the formula for l > 1:")
        print(calc_str)

    print("-" * 20)
    return result

def main():
    # Final expressions for the answer
    # d(n) is the divisor function.
    # delta_{i,j} is the Kronecker delta.
    # l = p_1^e_1 * ... * p_s^e_s
    ans_A = "d(l^2) - 1 + delta_{l,1}"
    ans_B = "((2*e_1+1)*(2*e_2+1)*...*(2*e_s+1)) - 1 + delta_{s,0}"
    
    # Print the answer in the requested format
    print(f"A)[{ans_A}] B)[{ans_B}]")
    
    # Example calculation for l=12
    result_12 = solve_and_explain(12)
    
    # Example calculation for l=1
    solve_and_explain(1)

    # Return the answer for l=12 in the final format
    print(f"\nThe example calculation for l=12 gives:")
    print(f"<<<{result_12}>>>")

if __name__ == "__main__":
    main()
