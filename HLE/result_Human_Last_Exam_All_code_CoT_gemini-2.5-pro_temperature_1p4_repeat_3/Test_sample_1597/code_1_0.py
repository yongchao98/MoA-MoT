import math

def calculate_allowed_pairs(n):
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a,b <= n.
    """
    # Step 1: Compute Mobius function mu up to n using a linear sieve.
    mu = [0] * (n + 1)
    lp = [0] * (n + 1)
    primes = []
    mu[1] = 1
    for i in range(2, n + 1):
        if lp[i] == 0:
            lp[i] = i
            mu[i] = -1
            primes.append(i)
        for p in primes:
            if p > lp[i] or i * p > n:
                break
            lp[i * p] = p
            if p == lp[i]:
                mu[i * p] = 0
            else:
                mu[i * p] = -mu[i]

    # Step 2: Calculate C(n), the number of coprime pairs (a,b) with 1 <= a,b <= n.
    coprime_pairs_count = 0
    for d in range(1, n + 1):
        term = mu[d] * (n // d)**2
        coprime_pairs_count += term

    # Step 3: Calculate the number of disallowed pairs.
    # Disallowed pairs are (a,b) where a>1, b>1, and gcd(a,b)=1.
    # This is C(n) minus the coprime pairs where a=1 or b=1.
    coprime_pairs_with_one = 2 * n - 1
    disallowed_pairs_count = coprime_pairs_count - coprime_pairs_with_one

    # Step 4: Calculate the total number of allowed pairs.
    total_pairs = n * n
    allowed_pairs_count = total_pairs - disallowed_pairs_count

    print(f"The calculation for N = {n}:")
    print(f"Total number of ordered pairs (a,b) is {n}x{n} = {total_pairs}.")
    print("A pair (a,b) is disallowed if a>1, b>1, and gcd(a,b)=1.")
    print(f"The number of coprime pairs (a,b) in [1,{n}]x[1,{n}] is C({n}) = {coprime_pairs_count}.")
    print(f"The number of coprime pairs with a=1 or b=1 is 2*{n}-1 = {coprime_pairs_with_one}.")
    print(f"Number of disallowed pairs = C({n}) - (2*{n}-1) = {coprime_pairs_count} - {coprime_pairs_with_one} = {disallowed_pairs_count}.")
    print(f"The number of allowed pairs = Total pairs - Disallowed pairs")
    print(f"= {total_pairs} - {disallowed_pairs_count} = {allowed_pairs_count}")

    return allowed_pairs_count

# Main execution
n = 1000
final_answer = calculate_allowed_pairs(n)