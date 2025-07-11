import math

def get_prime_factorization(num):
    """Returns a dictionary of the prime factors of num and their powers."""
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
        # If any prime factor is repeated (i.e., n is not square-free), mu(n) is 0.
        if factors[p] > 1:
            return 0
    
    # If n is square-free, mu(n) is (-1)^k where k is the number of distinct prime factors.
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

# --- Main Calculation ---
n = 10  # The order of the BCH expansion
k = 2   # The number of generators (X and Y)

print(f"Calculating the number of nonzero coefficients of order {n} using Witt's formula.")
print(f"L_n(k) = (1/n) * sum over d|n of (mu(d) * k^(n/d))")
print("-" * 40)

divisors = get_divisors(n)
print(f"For n = {n}, the divisors 'd' are: {', '.join(map(str, divisors))}\n")

total_sum = 0
calculation_terms = []

print("Calculating each term of the sum:")
for d in divisors:
    mu_d = mobius(d)
    exponent = n // d
    term_val = mu_d * (k**exponent)
    
    print(f"d = {d:2}:  mu({d}) * {k}^({n}/{d}) = {mu_d:2} * {k}^{exponent:<2} = {term_val:5}")
    
    total_sum += term_val
    calculation_terms.append(str(term_val))

sum_expression = " + ".join(f"({t})" if int(t) < 0 else t for t in calculation_terms)

print("\nFinal calculation:")
print(f"L_{n}({k}) = (1/{n}) * ({sum_expression})")
print(f"L_{n}({k}) = (1/{n}) * {total_sum}")

result = total_sum // n
print(f"L_{n}({k}) = {result}")

print(f"\nThe number of nonzero coefficients of order {n} is {result}.")
<<<99>>>