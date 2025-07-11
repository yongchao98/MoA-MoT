import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Computes the Mobius function mu(n).
    - mu(1) = 1
    - mu(n) = 0 if n has a squared prime factor.
    - mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    return (-1)**len(factors)

def get_divisors(n):
    """
    Returns a sorted list of all divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def calculate_bch_nonzero_coefficients():
    """
    Calculates and prints the number of nonzero coefficients of order n in the
    BCH expansion for k generators, following the problem's request.
    """
    n = 10
    k = 2

    print(f"To find the number of nonzero coefficients of order {n}, we use Witt's formula:")
    print(f"L_{n}({k}) = (1/{n}) * sum over d|{n} [ mu(d) * k^(n/d) ]\n")

    divisors = get_divisors(n)
    print(f"The divisors 'd' of {n} are: {divisors}\n")

    term_values = []
    
    # First, build and print the full formula with all numbers
    full_formula_parts = []
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_values.append(mu_d * (k**power))
        full_formula_parts.append(f"({mu_d} * {k}^{power})")
    
    print("Step 1: Substitute the values into the formula.")
    print(f"L_{n}({k}) = (1/{n}) * [ {' + '.join(full_formula_parts)} ]")

    # Second, print the formula with powers evaluated
    evaluated_powers_parts = []
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        evaluated_powers_parts.append(f"({mu_d} * {k**power})")
    
    print("\nStep 2: Evaluate the powers.")
    print(f"L_{n}({k}) = (1/{n}) * [ {' + '.join(evaluated_powers_parts)} ]")
    
    # Third, print the formula with multiplications performed
    multiplied_parts_str = " + ".join(map(str, term_values)).replace('+ -', '- ')
    
    print("\nStep 3: Perform the multiplications inside the brackets.")
    print(f"L_{n}({k}) = (1/{n}) * [ {multiplied_parts_str} ]")

    # Fourth, print the formula with the sum calculated
    total_sum = sum(term_values)
    print("\nStep 4: Calculate the sum inside the brackets.")
    print(f"L_{n}({k}) = (1/{n}) * [ {total_sum} ]")
    
    # Finally, print the result
    final_result = total_sum // n
    print("\nStep 5: Perform the final division to get the result.")
    print(f"L_{n}({k}) = {final_result}")

if __name__ == "__main__":
    calculate_bch_nonzero_coefficients()
