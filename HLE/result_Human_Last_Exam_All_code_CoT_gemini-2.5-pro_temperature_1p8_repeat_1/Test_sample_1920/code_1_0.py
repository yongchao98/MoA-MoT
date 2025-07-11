import math
from collections import defaultdict

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def get_divisors(n):
    """Returns a list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def count_chars_order_divides_k(p, k):
    """
    Counts primitive characters mod p whose order divides k.
    For a prime p, a character is primitive if its order > 1.
    """
    limit = math.gcd(k, p - 1)
    divs = get_divisors(limit)
    count = 0
    for j in divs:
        if j > 1: # Primitive characters mod p have order > 1
            count += phi(j)
    return count

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor d and order g.
    """
    d = 53599
    g = 6

    d_factors = get_prime_factorization(d)
    if any(exp > 1 for exp in d_factors.values()):
        print(f"Error: Conductor d={d} is not square-free. This method is for square-free conductors.")
        return

    g_prime_factors = list(get_prime_factorization(g).keys())
    
    total_count = 0
    
    explanation_lines = []

    # Inclusion-Exclusion principle on prime factors of g
    # Iterate through all subsets of the prime factors of g
    num_g_factors = len(g_prime_factors)
    for i in range(1 << num_g_factors):
        subset_product = 1
        num_in_subset = 0
        for j in range(num_g_factors):
            if (i >> j) & 1:
                subset_product *= g_prime_factors[j]
                num_in_subset += 1
        
        k = g // subset_product
        sign = (-1) ** num_in_subset
        
        # Calculate N_div(k), the number of tuples where lcm of orders divides k
        term_count = 1
        # Get count for a single p to build explanation string
        p_first = list(d_factors.keys())[0]
        count_for_one_p = count_chars_order_divides_k(p_first, k)
        
        # Calculate term product
        for p in d_factors:
            term_count *= count_chars_order_divides_k(p, k)
            
        total_count += sign * term_count

        # Build explanation
        op = "+" if sign == 1 else "-"
        # Since counts for each p are identical in this problem, we can use powers
        num_d_factors = len(d_factors)
        
        formula_part = f"({count_for_one_p})^{num_d_factors}"
        value_part = f"{term_count}"
        if i == 0: # First term has no preceding operator
             explanation_lines.append(f"N({k}) = {formula_part} = {value_part}")
        else:
             explanation_lines.append(f"{op} N({k}) = {op} {formula_part} = {op} {value_part}")


    print(f"Let d = {d} and g = {g}.")
    d_fact_str = ' * '.join(map(str, sorted(d_factors.keys())))
    print(f"The conductor d is square-free: d = {d_fact_str}.")
    print("The number of primitive characters of order g=6 is given by the inclusion-exclusion principle.")
    print("Let N(k) be the number of character tuples where the lcm of component orders divides k.")
    
    # We want to format the output equation clearly.
    # Total = N(6) - N(3) - N(2) + N(1)
    
    # For this specific problem, let's rearrange for clarity
    N_6_term = [s for s in explanation_lines if "N(6)" in s][0]
    N_3_term = [s for s in explanation_lines if "N(3)" in s][0]
    N_2_term = [s for s in explanation_lines if "N(2)" in s][0]
    N_1_term = [s for s in explanation_lines if "N(1)" in s][0]
    
    print("\nCalculation steps:")
    print(N_6_term)
    print(N_3_term.replace("- ", "")) # remove sign for separate line
    print(N_2_term.replace("- ", "")) # remove sign for separate line
    print(N_1_term.replace("+ ", "")) # remove sign for separate line

    final_eq = " ".join([N_6_term.split('=')[2].strip()] + [s.split('=')[1].strip() for s in [N_3_term, N_2_term, N_1_term]])
    print(f"\nFinal calculation: {final_eq} = {total_count}")

    print(f"\nThe number of primitive Dirichlet characters is {total_count}.")


solve()