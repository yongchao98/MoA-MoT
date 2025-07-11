def count_allowed_pairs(limit):
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= limit.
    An ordered pair (a,b) is allowed if for every primitive function f, f(ab) >= f(a)f(b).
    This condition holds if and only if a=1, b=1, or gcd(a,b) > 1.
    """

    # Step 1: Compute Mobius function mu(n) for n from 1 to limit using a sieve.
    mu = [0] * (limit + 1)
    lp = [0] * (limit + 1)
    primes = []
    mu[1] = 1

    for i in range(2, limit + 1):
        if lp[i] == 0:
            lp[i] = i
            primes.append(i)
            mu[i] = -1  # i is prime
        for p in primes:
            if p > lp[i] or i * p > limit:
                break
            lp[i * p] = p
            if p == lp[i]:  # p divides i/p
                mu[i * p] = 0
                break
            else: # p does not divide i/p
                mu[i * p] = -mu[i]

    # Step 2: Calculate the number of coprime pairs (a,b) with 1 <= a,b <= limit.
    # Formula: sum_{d=1 to limit} mu(d) * floor(limit/d)^2
    coprime_pairs_count = 0
    for d in range(1, limit + 1):
        term = mu[d] * (limit // d)**2
        coprime_pairs_count += term

    # Step 3: Calculate the number of not-allowed pairs.
    # Not-allowed pairs are those with a > 1, b > 1, and gcd(a,b) = 1.
    # The coprime_pairs_count includes pairs where a=1 or b=1.
    # Number of coprime pairs with a=1 or b=1 is (2 * limit - 1).
    # (1,1), (1,2)...(1,limit) -> limit pairs
    # (2,1), (3,1)...(limit,1) -> limit-1 pairs
    not_allowed_pairs_count = coprime_pairs_count - (2 * limit - 1)

    # Step 4: Calculate the total number of allowed pairs.
    total_pairs = limit * limit
    allowed_pairs_count = total_pairs - not_allowed_pairs_count
    
    print(f"The number of allowed pairs (a,b) with 1 <= a,b <= 1000 is calculated as follows:")
    print(f"Total number of pairs = {limit} * {limit} = {total_pairs}")
    print(f"Number of coprime pairs (1 <= a,b <= {limit}) = {coprime_pairs_count}")
    print(f"Number of not-allowed pairs (a>1, b>1, gcd(a,b)=1) = {coprime_pairs_count} - (2*{limit}-1) = {not_allowed_pairs_count}")
    print(f"Number of allowed pairs = {total_pairs} - {not_allowed_pairs_count} = {allowed_pairs_count}")

if __name__ == '__main__':
    count_allowed_pairs(1000)