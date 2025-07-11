def get_prime_factorization(num):
    """
    Returns the prime factorization of a number as a dictionary.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while temp_num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Möbius function μ(n).
    μ(1) = 1
    μ(n) = 0 if n has a squared prime factor.
    μ(n) = (-1)^k if n is a product of k distinct primes.
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    return (-1)**len(factors)

def solve_bch_problem():
    """
    Calculates the number of nonzero coefficients of a given order in the
    Baker-Campbell-Hausdorff expansion and prints the detailed steps.
    """
    n = 10
    k = 2

    print("The number of nonzero coefficients of order 10 is calculated using Witt's formula.")
    print("This formula gives the dimension of the degree-n component of the free Lie algebra on k generators.\n")
    print(f"Witt's Formula: L_k(n) = (1/n) * Σ_{{d|n}} μ(d) * k^(n/d)\n")
    
    # Find divisors of n
    divisors = sorted([d for d in range(1, n + 1) if n % d == 0])
    
    print(f"For our case, n = {n}, k = {k}.")
    print(f"The divisors 'd' of {n} are: {divisors}\n")
    
    total_sum = 0
    
    # Build the equation parts for printing
    symbolic_terms = []
    numeric_terms_expanded = []
    numeric_terms_simple = []
    
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_value = mu_d * (k ** power)
        total_sum += term_value
        
        symbolic_terms.append(f"μ({d})*({k}^({n}/{d}))")
        numeric_terms_expanded.append(f"({mu_d})*({k**power})")
        numeric_terms_simple.append(str(term_value))

    # Print step 1: The symbolic formula
    print("Step 1: Substitute n, k, and the divisors into the formula.")
    equation_step1 = f"L_{k}({n}) = (1/{n}) * (" + " + ".join(symbolic_terms) + ")"
    print(equation_step1)
    
    # Print step 2: The formula with calculated values for μ(d) and k^(n/d)
    print("\nStep 2: Calculate μ(d) and k^(n/d) for each divisor d.")
    equation_step2 = f"L_{k}({n}) = (1/{n}) * (" + " + ".join(numeric_terms_expanded) + ")"
    print(equation_step2.replace('+ (-', '- (')) # for cleaner display

    # Print step 3: The formula with evaluated terms
    print("\nStep 3: Evaluate each term in the sum.")
    simple_sum_str = " ".join(f"- {s[1:]}" if s.startswith('-') else f"+ {s}" for s in numeric_terms_simple).lstrip('+ ')
    equation_step3 = f"L_{k}({n}) = (1/{n}) * ({simple_sum_str})"
    print(equation_step3)
    
    # Print step 4: The sum of the terms
    print("\nStep 4: Calculate the total sum inside the parenthesis.")
    equation_step4 = f"L_{k}({n}) = (1/{n}) * ({total_sum})"
    print(equation_step4)
    
    # Print step 5: The final result
    final_result = total_sum // n
    print("\nStep 5: Perform the final division.")
    equation_step5 = f"L_{k}({n}) = {final_result}"
    print(equation_step5)
    
    print(f"\nThe number of nonzero coefficients of order 10 is {final_result}.")

solve_bch_problem()
<<<99>>>