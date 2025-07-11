import collections

def get_prime_factorization(n):
    """
    Returns a dictionary where keys are prime factors and values are their exponents.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling 100 six-sided dice.
    """
    num_dice = 100
    num_sides = 6

    print(f"Problem: Find the largest number of mutually independent events for rolling {num_dice} {num_sides}-sided dice.")
    print("-" * 80)
    print("The solution is based on the prime factorization of the total number of outcomes in the sample space.")
    
    # Step 1: Find prime factorization of the number of sides.
    side_factors = get_prime_factorization(num_sides)
    
    factor_str_parts = []
    for p, e in side_factors.items():
        factor_str_parts.append(f"{p}^{e}" if e > 1 else str(p))
    print(f"\n1. The prime factorization of the number of sides ({num_sides}) is: {' * '.join(factor_str_parts)}")

    # Step 2: Determine the prime factorization of the total number of outcomes.
    print(f"\n2. The total number of outcomes is {num_sides}^{num_dice}, so its prime factorization is:")
    
    n_factor_str_parts = []
    exponent_terms = []
    total_m = 0
    
    for p, e in sorted(side_factors.items()):
        total_exponent = e * num_dice
        n_factor_str_parts.append(f"{p}^{total_exponent}")
        exponent_terms.append(str(total_exponent))
        total_m += total_exponent
        
    print(f"   N = ({' * '.join(factor_str_parts)})^{num_dice} = {' * '.join(n_factor_str_parts)}")
    
    # Step 3: Sum the exponents to find the maximum number of independent events (m).
    print("\n3. The maximum number of mutually independent events (m) is the sum of the exponents in the prime factorization of N.")
    
    # Output the final equation with each number
    final_equation = f"   m = {' + '.join(exponent_terms)}"
    print(final_equation)
    
    final_result = f"   m = {total_m}"
    print(final_result)
    
    print("-" * 80)


solve_dice_problem()
<<<200>>>