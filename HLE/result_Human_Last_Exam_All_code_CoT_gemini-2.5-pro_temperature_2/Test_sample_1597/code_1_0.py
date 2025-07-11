import math

def count_allowed_pairs(N):
    """
    Calculates the number of 'allowed' pairs (a,b) with 1 <= a,b <= N.

    An ordered pair (a,b) is allowed if a=1, b=1, or gcd(a,b) > 1.
    We count the number of disallowed pairs (a>1, b>1, gcd(a,b)=1) and subtract from the total.
    """

    # Step 1: Compute Mobius function values up to N using a sieve
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
            else:
                mu[i * p] = -mu[i]

    # Step 2: Count all pairs (a,b) with 1 <= a,b <= N and gcd(a,b) = 1
    total_coprime_pairs = 0
    for d in range(1, N + 1):
        term = mu[d] * (N // d) * (N // d)
        total_coprime_pairs += term
        
    # Step 3: Count coprime pairs where a=1 or b=1
    # This is (1,k) for 1<=k<=N, and (k,1) for 1<=k<=N.
    # The pair (1,1) is counted in both sets.
    # Since gcd(1,k)=1 for any k, all these pairs are coprime.
    coprime_pairs_with_one = 2 * N - 1

    # Step 4: Calculate the number of disallowed pairs
    # These are pairs (a,b) where a>1, b>1, and gcd(a,b)=1
    disallowed_count = total_coprime_pairs - coprime_pairs_with_one

    # Step 5: Calculate the number of allowed pairs
    total_pairs = N * N
    allowed_count = total_pairs - disallowed_count

    print(f"For N = {N}:")
    print(f"Total number of pairs = {N} * {N} = {total_pairs}")
    print(f"Total number of coprime pairs (1 <= a,b <= {N}) = {total_coprime_pairs}")
    print(f"Number of disallowed pairs (a>1, b>1, gcd(a,b)=1) = {total_coprime_pairs} - {coprime_pairs_with_one} = {disallowed_count}")
    print(f"Total number of allowed pairs = {total_pairs} - {disallowed_count} = {allowed_count}")

# Set N = 1000 for the problem
N = 1000
count_allowed_pairs(N)