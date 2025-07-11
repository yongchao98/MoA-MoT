import math

def get_prime_factorization(n):
    """
    Computes the prime factorization of a positive integer n.
    Returns a dictionary where keys are prime factors and values are their exponents.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve_for_l(l):
    """
    Calculates and prints the cardinalities of U_l and T_l for a given l.
    """
    if not isinstance(l, int) or l <= 0:
        print("Error: l must be a positive integer.")
        return

    # Part A
    factors = get_prime_factorization(l)
    s = len(factors)
    card_U_l = 2**s
    
    # Part B
    e = list(factors.values())
    if l == 1:
        card_T_l = 1
    else:
        prod = 1
        for val in e:
            prod *= (1 + 2 * val)
        card_T_l = prod - 1
        
    # --- Final Answer Output ---
    
    # Part A answer
    print("A) |U_l| has the following expression in terms of s (the number of distinct prime factors):")
    final_eq_A = f"2^s"
    print(final_eq_A)
    print(f"For l = {l}, s = {s}, so the value is {card_U_l}")
    print("-" * 20)
    
    # Part B answer
    print("B) |T_l| can be calculated as follows:")
    if l == 1:
        print("For l = 1, the value is 1.")
    else:
        print("The expression is: (product over i=1 to s of (1+2*e_i)) - 1")
        val_strs = [f"(1+2*{v})" for v in e]
        # This part of the code prints each number in the final equation as requested.
        final_eq_B = f"({' * '.join(val_strs)}) - 1"
        print(f"For l = {l}, the exponents are {e}.")
        print(f"The calculation is: {final_eq_B} = {card_T_l}")

# Example for l = 12
l_example = 12
solve_for_l(l_example)

# Example for l = 1
print("\n")
l_example = 1
solve_for_l(l_example)