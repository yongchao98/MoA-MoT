import math

def mobius_sieve(n):
    """
    Computes the Mobius function mu(k) for all k from 1 to n.
    mu(k) = 1 if k is a square-free positive integer with an even number of prime factors.
    mu(k) = -1 if k is a square-free positive integer with an odd number of prime factors.
    mu(k) = 0 if k has a squared prime factor.
    """
    mu = [0] * (n + 1)
    is_prime = [True] * (n + 1)
    primes = []
    mu[1] = 1
    
    for i in range(2, n + 1):
        if is_prime[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_prime[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]
    return mu

def solve():
    """
    Counts the number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    A pair is allowed if gcd(a,b) is square-free.
    """
    N = 1000
    limit = int(math.sqrt(N))
    mu = mobius_sieve(limit)

    disallowed_count = 0
    
    # Using the Principle of Inclusion-Exclusion to count disallowed pairs.
    # A pair (a,b) is disallowed if gcd(a,b) is divisible by p^2 for some prime p.
    # The number of such pairs is given by: Sum_{k>1, k is sq-free} -mu(k) * (floor(N/k^2))^2
    
    print("The number of disallowed pairs is the sum of the following terms (using inclusion-exclusion):")
    
    # The final equation is Disallowed = Sum -mu(k) * floor(N/k^2)^2
    # We print each non-zero term in this sum.
    equation_parts = []

    for k in range(2, limit + 1):
        if mu[k] == 0:
            continue
        
        term = -mu[k] * ((N // (k * k)) ** 2)
        disallowed_count += term
        
        if term > 0:
            equation_parts.append(f"+ {term}")
        else:
            equation_parts.append(f"- {abs(term)}")
    
    # Let's print the equation terms properly
    first_term = equation_parts[0].replace("+ ", "")
    print(f"{first_term} {' '.join(equation_parts[1:])}")
    
    total_pairs = N * N
    allowed_count = total_pairs - disallowed_count
    
    print(f"\nThe calculation is: {disallowed_count}")
    print(f"Total number of pairs from 1 to {N} is {total_pairs}.")
    print(f"Number of disallowed pairs where gcd(a,b) is not square-free is {disallowed_count}.")
    print(f"The number of allowed pairs is {total_pairs} - {disallowed_count} = {allowed_count}.")
    
solve()