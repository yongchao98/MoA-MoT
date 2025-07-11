import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        if temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        else:
            d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    # If any prime factor has a power greater than 1, mu(n) is 0
    if any(p > 1 for p in factors.values()):
        return 0
    
    # Otherwise, it's (-1)^k where k is the number of distinct prime factors
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def main():
    """
    Calculates the number of nonzero coefficients of order n in BCH expansion
    and prints the step-by-step formula.
    """
    # Order of the expansion
    n = 10
    # Number of generators (X and Y)
    k = 2

    # Step 1: Find divisors
    divisors = get_divisors(n)

    # Step 2 & 3: Calculate terms for Witt's formula
    mu_values = [mobius(d) for d in divisors]
    powers = [n // d for d in divisors]
    term_values = [m * (k**p) for m, p in zip(mu_values, powers)]
    
    # Step 4: Sum the terms
    total_sum = sum(term_values)
    
    # Step 5: Get final result
    result = total_sum // n

    # Step 6: Print the detailed calculation
    print(f"The number of nonzero coefficients of order {n} in the Baker-Campbell-Hausdorff expansion is given by Witt's formula, which evaluates to {result}.")
    print("The calculation is as follows:")

    # Build the first line of the equation: result = (1/n) * [ (mu(d1))*k^(n/d1) + ... ]
    part1_terms = []
    for mu, d, p in zip(mu_values, divisors, powers):
        part1_terms.append(f"({mu})*{k}^{p}")
    part1_str = " + ".join(part1_terms)
    print(f"{result} = (1/{n}) * [ {part1_str} ]")

    # Build the second line of the equation: result = (1/n) * [ term1 + term2 + ... ]
    part2_str = str(term_values[0])
    for term in term_values[1:]:
        if term >= 0:
            part2_str += f" + {term}"
        else:
            part2_str += f" - {-term}"
    print(f"{result} = (1/{n}) * [ {part2_str} ]")

    # Build the third line: result = (1/n) * [ sum ]
    print(f"{result} = (1/{n}) * [ {total_sum} ]")
    
    # Final answer
    print(f"{result} = {result}")

if __name__ == "__main__":
    main()