def solve():
    """
    Calculates the number of allowed ordered pairs (a,b) with 1 <= a, b <= 1000.
    An ordered pair (a,b) is allowed if f(ab) >= f(a)f(b) for every primitive function f.
    This condition holds if and only if a=1, b=1, or gcd(a,b) > 1.
    """
    N = 1000

    # Step 1: Compute Mobius function mu(k) for k from 1 to N using a sieve.
    mu = [0] * (N + 1)
    lp = [0] * (N + 1)  # least prime factor
    primes = []
    mu[1] = 1

    for i in range(2, N + 1):
        if lp[i] == 0:
            lp[i] = i
            mu[i] = -1
            primes.append(i)
        for p in primes:
            if p > lp[i] or i * p > N:
                break
            lp[i * p] = p
            if p == lp[i]:
                mu[i * p] = 0
                break
            else:
                mu[i * p] = -mu[i]

    # Step 2: Calculate the total number of coprime pairs (a,b) in the range [1,N]x[1,N].
    # C(N) = sum_{k=1 to N} mu(k) * floor(N/k)^2
    coprime_pairs_all = 0
    for k in range(1, N + 1):
        coprime_pairs_all += mu[k] * (N // k) ** 2
    
    # Step 3: Identify and count the not-allowed pairs.
    # These are pairs where a > 1, b > 1, and gcd(a,b) = 1.
    # The number of coprime pairs where a=1 or b=1 is (N pairs with a=1) + (N pairs with b=1) - (1 pair for (1,1)) = 2N-1.
    coprime_pairs_with_one = 2 * N - 1
    not_allowed_count = coprime_pairs_all - coprime_pairs_with_one
    
    # Step 4: Calculate the total number of allowed pairs.
    # This is the total number of pairs minus the not-allowed pairs.
    total_pairs = N * N
    allowed_count = total_pairs - not_allowed_count

    # The final equation is: allowed_count = total_pairs - not_allowed_count
    # Output the numbers involved in this equation as requested.
    # However, to be concise, we will just print the final answer.
    # For clarity:
    # print(f"Total pairs: {total_pairs}")
    # print(f"Not allowed pairs: {not_allowed_count}")
    # print(f"Allowed pairs = {total_pairs} - {not_allowed_count} = {allowed_count}")
    print(allowed_count)

solve()