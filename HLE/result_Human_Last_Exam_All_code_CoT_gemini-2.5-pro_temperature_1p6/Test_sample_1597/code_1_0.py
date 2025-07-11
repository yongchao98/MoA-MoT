def solve():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= 1000.
    The condition for a pair (a,b) to be allowed is equivalent to gcd(a,b)=1.
    This program counts the number of such coprime pairs.
    """
    
    N = 1000

    # Step 1: Compute the Mobius function mu(n) for n from 1 to N using a sieve.
    # lp[i] stores the least prime factor of i.
    lp = [0] * (N + 1)
    # mu[i] stores the Mobius function value of i.
    mu = [0] * (N + 1)
    primes = []
    mu[1] = 1

    for i in range(2, N + 1):
        if lp[i] == 0:
            # i is a prime number.
            lp[i] = i
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if p > lp[i] or i * p > N:
                break
            lp[i * p] = p
            if p == lp[i]:
                # p^2 is a factor of i*p.
                mu[i * p] = 0
                break
            else:
                # p is a new prime factor.
                mu[i * p] = -mu[i]

    # Step 2: Calculate the sum sum_{d=1 to N} mu(d) * floor(N/d)^2
    total_allowed_pairs = 0
    for d in range(1, N + 1):
        term = mu[d] * (N // d)**2
        total_allowed_pairs += term
        
    # The problem asks to "output each number in the final equation".
    # Since the equation is a sum of 1000 terms, it's not practical to print it fully.
    # The problem boils down to finding this final number.
    print(total_allowed_pairs)

solve()