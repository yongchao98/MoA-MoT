def count_allowed_pairs(N):
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= N.
    An ordered pair (a,b) is allowed if a=1, b=1, or gcd(a,b) > 1.
    """

    # Step 1: Compute Mobius function mu(d) for d from 1 to N using a sieve.
    mu = [0] * (N + 1)
    lp = [0] * (N + 1)
    primes = []
    mu[1] = 1

    for i in range(2, N + 1):
        if lp[i] == 0:
            lp[i] = i
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if p > lp[i] or i * p > N:
                break
            lp[i * p] = p
            if p == lp[i]:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]

    # Step 2: Calculate C(N), the number of coprime pairs (a,b) with 1 <= a,b <= N.
    # C(N) = sum_{d=1 to N} mu(d) * floor(N/d)^2
    num_coprime_pairs = 0
    for d in range(1, N + 1):
        num_coprime_pairs += mu[d] * (N // d)**2

    # Step 3: Calculate the number of disallowed pairs.
    # These are pairs (a,b) where a>1, b>1, and gcd(a,b)=1.
    # This is C(N) minus the coprime pairs involving 1.
    # The coprime pairs involving 1 are (1, 1..N) and (2..N, 1), totaling 2N-1 pairs.
    num_disallowed_pairs = num_coprime_pairs - (2 * N - 1)

    # Step 4: Calculate the number of allowed pairs.
    # This is the total number of pairs minus the disallowed ones.
    total_pairs = N * N
    num_allowed_pairs = total_pairs - num_disallowed_pairs
    
    print(f"The number of ordered pairs of positive integers (a,b) such that 1 <= a,b <= 1000 that are allowed is: {num_allowed_pairs}")

# Run the calculation for N=1000
count_allowed_pairs(1000)