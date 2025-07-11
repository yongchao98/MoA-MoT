import math

def count_allowed_pairs(N):
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a, b <= N.
    """
    # Step 1: Find all primes up to N using a sieve
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(N)) + 1):
        if is_prime[i]:
            for multiple in range(i * i, N + 1, i):
                is_prime[multiple] = False
    primes = [i for i, p_bool in enumerate(is_prime) if p_bool]

    # Step 2: For each prime p, calculate k_p = floor(log_p(N))
    # This is the number of powers of p (p^1, p^2, ...) that are <= N.
    k_p_values = []
    for p in primes:
        if p == 0: continue # Should not happen with our sieve
        k_p = int(math.log(N) / math.log(p))
        if k_p > 0:
            k_p_values.append(k_p)

    # Step 3: Calculate sum_kp and sum_kp_sq
    # sum_kp is the total number of prime powers <= N (excluding 1)
    # The number of disallowed pairs = (sum k_p)^2 - (sum k_p^2)
    sum_kp = sum(k_p_values)
    sum_kp_sq = sum(k * k for k in k_p_values)
    
    # A pair (a,b) is disallowed if a=p^k and b=q^l for p != q.
    disallowed_pairs_count = sum_kp * sum_kp - sum_kp_sq
    
    # Step 4: Calculate total pairs and allowed pairs
    total_pairs_count = N * N
    allowed_pairs_count = total_pairs_count - disallowed_pairs_count
    
    print(f"Total number of pairs (a,b) where 1 <= a,b <= {N} is: {N} * {N} = {total_pairs_count}")
    print(f"Let k_p be the number of powers of a prime p which are less than or equal to {N}.")
    print(f"The sum of all k_p is Sum(k_p) = {sum_kp}")
    print(f"The sum of all k_p^2 is Sum(k_p^2) = {sum_kp_sq}")
    print(f"The number of disallowed pairs is (Sum(k_p))^2 - Sum(k_p^2) = {sum_kp}^2 - {sum_kp_sq} = {disallowed_pairs_count}")
    print(f"The number of allowed pairs is Total_Pairs - Disallowed_Pairs = {total_pairs_count} - {disallowed_pairs_count} = {allowed_pairs_count}")
    
    # Final answer in required format.
    # print(f"\n<<<{allowed_pairs_count}>>>")

if __name__ == '__main__':
    N = 1000
    count_allowed_pairs(N)
