def solve():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    The condition for a pair (a,b) to be allowed is gcd(a,b) = 1.
    The number of such pairs is calculated using the formula:
    Sum_{d=1 to N} [mu(d) * floor(N/d)^2], where N=1000.
    """
    N = 1000

    # Step 1: Compute Mobius function mu(d) for d from 1 to N using a sieve.
    mu = [0] * (N + 1)
    is_prime = [True] * (N + 1)
    primes = []
    mu[1] = 1

    for i in range(2, N + 1):
        if is_prime[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > N:
                break
            is_prime[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]

    # Step 2: Calculate the sum.
    total_allowed_pairs = 0
    for d in range(1, N + 1):
        term = mu[d] * (N // d) * (N // d)
        total_allowed_pairs += term

    print(total_allowed_pairs)

solve()