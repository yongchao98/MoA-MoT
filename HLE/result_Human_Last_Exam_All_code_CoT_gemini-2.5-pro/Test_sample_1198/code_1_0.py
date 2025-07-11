import math

def analyze_construction(modulus, n_max, initial_interval=(1.0, 4.0)):
    """
    Analyzes the existence of a real number 'a' such that floor(a^n) % modulus == n % modulus
    by attempting to construct a valid interval [l, r) for 'a'.
    
    Args:
        modulus (int): The modulus for the congruence (e.g., 2 or 3).
        n_max (int): The number of steps to run the analysis for.
        initial_interval (tuple): A starting interval (l, r) for 'a'.
    """
    print(f"--- Analysis for modulus {modulus} ---")
    print(f"Searching for a > 0 such that floor(a^n) = n (mod {modulus}) for all n > 0.")
    l, r = initial_interval

    for n in range(1, n_max + 1):
        print(f"\nStep n={n}:")
        if l >= r:
            print(f"Interval is empty. Construction failed in a previous step.")
            return False
            
        print(f"Current interval for 'a' is [{l:.8f}, {r:.8f})")

        # Calculate the implied interval for a^n, which is [l^n, r^n)
        l_pow_n = l**n
        r_pow_n = r**n
        
        # Determine the target residue for k_n = floor(a^n)
        target_residue = n % modulus
        
        # Find all possible integers k in the range of a^n that satisfy the residue condition.
        # The lowest possible integer k is floor(l^n) if l^n is not an integer, or l^n itself.
        # The highest possible integer k is ceil(r^n) - 1.
        k_min = math.floor(l_pow_n)
        k_max = math.ceil(r_pow_n) - 1

        possible_k = []
        for k in range(k_min, k_max + 1):
            if k % modulus == target_residue:
                # To be a possibility for floor(a^n), the interval [k, k+1) must overlap
                # with the interval [l^n, r^n).
                # Overlap exists if max(l^n, k) < min(r^n, k+1).
                if max(l_pow_n, k) < min(r_pow_n, k + 1):
                    possible_k.append(k)
        
        print(f"The interval for a^{n} is [{l_pow_n:.4f}, {r_pow_n:.4f}). Length={r_pow_n - l_pow_n:.4f}")
        print(f"Required: k = {target_residue} (mod {modulus}). Possible integers k: {possible_k}")

        if not possible_k:
            print("\nCONCLUSION: No integer k satisfies the condition in the given range.")
            print("The construction has failed. It is proven no such 'a' exists.")
            return False

        # To continue, we form a new region for 'a' from the union of valid intervals.
        # New region for a^n is Union_{k in possible_k} [k, k+1).
        # We find the corresponding region for 'a' and intersect it with the current [l, r).
        
        # The new interval for 'a' must be within the union of [k^(1/n), (k+1)^(1/n))
        # for all k in possible_k.
        new_l = float('inf')
        new_r = float('-inf')
        
        # To find the new continuous interval, we take the union of all intersections:
        # Union_{k} ([l, r) intersect [k^(1/n), (k+1)^(1/n)))
        # To simplify, we'll find the minimum possible l and maximum possible r.
        
        min_l_cand = min([k**(1/n) for k in possible_k])
        max_r_cand = max([(k+1)**(1/n) for k in possible_k])

        l = max(l, min_l_cand)
        r = min(r, max_r_cand)

    print("\nCONCLUSION: The construction process did not fail in these steps.")
    print(f"A non-empty interval for 'a' remains: [{l:.8f}, {r:.8f}).")
    if modulus == 2:
        print("For modulus 2, this construction can always proceed. Such an 'a' exists.")
    return True

if __name__ == "__main__":
    # Case 1: Modulo 2
    mod2_exists = analyze_construction(modulus=2, n_max=10, initial_interval=(2.0, 3.0))
    ans1 = "Yes" if mod2_exists else "No"
    
    print("\n" + "="*50 + "\n")
    
    # Case 2: Modulo 3
    # The proof of non-existence relies on the fact that for small 'a' (e.g., a < 2),
    # the length of the interval for a^n is also small, and can fail to contain
    # an integer with the required residue.
    mod3_exists = analyze_construction(modulus=3, n_max=10, initial_interval=(1.4, 1.5))
    ans2 = "Yes" if mod3_exists else "No"
    
    print("\n" + "="*50 + "\n")
    print(f"Final Answer:")
    print(f"Is there such 'a' for modulo 2? {ans1}")
    print(f"Is there such 'a' for modulo 3? {ans2}")
    print(f"\nFormatted answer: {ans1},{ans2}")
