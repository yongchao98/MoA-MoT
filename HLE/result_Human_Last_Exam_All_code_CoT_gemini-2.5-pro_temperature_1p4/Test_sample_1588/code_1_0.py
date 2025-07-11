import math

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary."""
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
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        # If any prime factor is squared, mu(n) = 0
        if factors[p] > 1:
            return 0
    
    # Otherwise, mu(n) = (-1)^k where k is the number of distinct prime factors
    return (-1) ** len(factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of a number."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_problem():
    """
    Calculates the number of nonzero coefficients of a given order in the
    Baker-Campbell-Hausdorff expansion using Witt's formula and prints the steps.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators (X and Y)

    print(f"The number of nonzero coefficients of order {n} is given by Witt's formula:")
    print(f"L_k(n) = (1/n) * Sum_{{d|n}} mu(d) * k^(n/d)\n")
    print(f"For our case, n={n} and k={k}:")
    
    divs = get_divisors(n)
    mu_values = [mobius(d) for d in divs]
    
    # 1. Print the formula with symbols
    sum_str = " + ".join([f"mu({d})*2^({n}/{d})" for d in divs])
    print(f"L_2({n}) = (1/{n}) * ({sum_str})")

    # 2. Print the formula with mu values substituted
    mu_val_str_parts = [f"({mu})*2^({n//d})" for mu, d in zip(mu_values, divs)]
    mu_val_str = " + ".join(mu_val_str_parts)
    print(f"L_2({n}) = (1/{n}) * ({mu_val_str})")

    # 3. Print the formula with powers calculated
    pow_val_str_parts = [f"({mu})*{k**(n//d)}" for mu, d in zip(mu_values, divs)]
    pow_val_str = " + ".join(pow_val_str_parts)
    print(f"L_2({n}) = (1/{n}) * ({pow_val_str})")

    # 4. Print the formula with each term calculated
    terms = [mu * (k**(n//d)) for mu, d in zip(mu_values, divs)]
    term_str = " + ".join([str(t) if t >= 0 else f"({t})" for t in terms])
    print(f"L_2({n}) = (1/{n}) * ({term_str})")

    # 5. Print the formula with the sum calculated
    total_sum = sum(terms)
    print(f"L_2({n}) = (1/{n}) * ({total_sum})")
    
    # 6. Print the final result
    result = total_sum // n
    print(f"L_2({n}) = {result}\n")
    
    print(f"The number of nonzero coefficients of order 10 is {result}.")

if __name__ == '__main__':
    solve_bch_problem()