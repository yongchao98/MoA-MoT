import math

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
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

def mobius_mu(n):
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0  # n has a squared prime factor
            
    # n is square-free
    return (-1)**len(factors)

def main():
    """
    Calculates the number of nonzero BCH coefficients of order n for k generators.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    print(f"The number of nonzero coefficients of order n in the Baker-Campbell-Hausdorff expansion")
    print(f"for k generators is given by Witt's formula:")
    print(f"N(n) = (1/n) * sum_{{d|n}} mu(d) * k^(n/d)\n")
    
    print(f"For this problem, n = {n} and k = {k}.")
    
    divisors = get_divisors(n)
    print(f"The positive divisors 'd' of {n} are: {', '.join(map(str, divisors))}\n")

    print("Calculating each term 'mu(d) * k^(n/d)' of the sum:")
    
    total_sum = 0
    term_expressions = []
    for d in divisors:
        mu_d = mobius_mu(d)
        power = n // d
        term_val = mu_d * (k**power)
        total_sum += term_val
        
        term_str = f"mu({d}) * {k}^({n}/{d})"
        value_str = f"{mu_d} * {k}^{power}"
        print(f"d={d}:  {term_str:<18} = {value_str:<12} = {term_val}")
        term_expressions.append(f"{term_val}")
        
    print("\nNow, substituting these values into the formula:")
    sum_str = " + ".join(term_expressions).replace("+ -", "- ")
    
    print(f"N({n}) = (1/{n}) * ({sum_str})")
    print(f"N({n}) = (1/{n}) * ({total_sum})")
    
    result = total_sum // n
    print(f"N({n}) = {result}")

if __name__ == "__main__":
    main()