import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of a number.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    """
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    # Check for squared prime factors
    for p in factors:
        if factors[p] > 1:
            return 0
    
    # If square-free, return (-1)^k where k is the number of distinct prime factors
    return (-1)**len(factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of a number.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_order_count(n, k):
    """
    Calculates the number of nonzero coefficients for order n and k generators
    and prints the step-by-step calculation.
    """
    print(f"Calculating the number of nonzero coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: L_k(n) = (1/n) * sum_{d|n} mu(d) * k^(n/d)")
    
    divisors = get_divisors(n)
    print(f"\nThe positive divisors 'd' of n={n} are: {divisors}")
    
    sum_terms_values = []
    sum_terms_str = []
    
    for d in divisors:
        mu_d = mobius(d)
        power_val = k**(n // d)
        term_value = mu_d * power_val
        
        sum_terms_values.append(term_value)
        sum_terms_str.append(f"({mu_d} * {k}^{n//d})")
        
    total_sum = sum(sum_terms_values)
    result = total_sum // n
    
    # Build the equation string
    # e.g., (1/10) * [(1 * 1024) + (-1 * 32) + (-1 * 4) + (1 * 2)]
    sum_calc_str = " + ".join(sum_terms_str)
    
    print(f"\nStep 1: Construct the sum for L_{k}({n})")
    print(f"L_{k}({n}) = (1/{n}) * [{sum_calc_str}]")

    # Build the next step of calculation string
    # e.g., (1/10) * [1024 - 32 - 4 + 2]
    values_str_list = [str(v) for v in sum_terms_values]
    calc_str = " + ".join(values_str_list).replace("+ -", "- ")
    
    print(f"\nStep 2: Evaluate each term in the sum")
    print(f"L_{k}({n}) = (1/{n}) * [{calc_str}]")
    
    print(f"\nStep 3: Calculate the total sum")
    print(f"L_{k}({n}) = (1/{n}) * {total_sum}")

    print(f"\nStep 4: Perform the final division")
    print(f"L_{k}({n}) = {total_sum} / {n} = {result}")

    print(f"\nThus, the number of nonzero coefficients of order {n} is {result}.")
    
    
# Parameters for the problem
n_order = 10
k_generators = 2

solve_bch_order_count(n_order, k_generators)