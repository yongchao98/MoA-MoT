import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num."""
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
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coeffs(n, k):
    """
    Calculates the number of nonzero coefficients of order n in the BCH expansion
    for k generators using Witt's formula and prints the calculation steps.
    """
    print(f"Calculating the number of nonzero coefficients of order {n} for {k} generators.")
    print("This is given by Witt's formula: L_n(k) = (1/n) * sum_{d|n} (μ(d) * k^(n/d))")
    print("-" * 20)

    divs = get_divisors(n)
    print(f"Step 1: Find the divisors of n = {n}")
    print(f"Divisors are: {divs}")
    print()

    mu_values = {d: mobius(d) for d in divs}
    print(f"Step 2: Calculate the Möbius function μ(d) for each divisor d")
    for d in divs:
        print(f"μ({d}) = {mu_values[d]}")
    print()

    print(f"Step 3: Apply Witt's formula for n={n} and k={k}")
    
    # Build the formula string
    formula_str = f"L_{n}({k}) = (1/{n}) * ("
    terms_str = []
    for d in divs:
        terms_str.append(f"μ({d})*{k}^({n}/{d})")
    formula_str += " + ".join(terms_str) + ")"
    print(formula_str)

    # Build the values string
    values_str = f"L_{n}({k}) = (1/{n}) * ("
    terms_vals = []
    for d in divs:
        mu_d = mu_values[d]
        # Use parenthesis for negative numbers
        mu_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        terms_vals.append(f"{mu_str}*{k}^{int(n/d)}")
    values_str += " + ".join(terms_vals) + ")"
    print(values_str)
    
    # Build the evaluated powers string
    eval_powers_str = f"L_{n}({k}) = (1/{n}) * ("
    terms_eval_powers = []
    total_sum = 0
    for d in divs:
        mu_d = mu_values[d]
        power_val = k**int(n/d)
        term = mu_d * power_val
        total_sum += term
        
        mu_str = f"({mu_d})" if mu_d < 0 else str(mu_d)
        terms_eval_powers.append(f"{mu_str}*{power_val}")
    eval_powers_str += " + ".join(terms_eval_powers) + ")"
    print(eval_powers_str)

    # Build the final calculation string
    final_calc_str = f"L_{n}({k}) = (1/{n}) * ("
    terms_final_calc = []
    for d in divs:
        mu_d = mu_values[d]
        power_val = k**int(n/d)
        term = mu_d * power_val
        terms_final_calc.append(str(term))
    final_calc_str += " + ".join(terms_final_calc).replace("+ -", "- ") + ")"
    print(final_calc_str)

    print(f"L_{n}({k}) = (1/{n}) * ({total_sum})")
    
    result = total_sum // n
    print(f"L_{n}({k}) = {result}")
    print("-" * 20)
    print(f"The number of nonzero coefficients of order {n} is {result}.")
    return result

# For the specific problem: n=10, k=2
n_order = 10
k_generators = 2
final_answer = solve_bch_coeffs(n_order, k_generators)
print(f"<<<{final_answer}>>>")
