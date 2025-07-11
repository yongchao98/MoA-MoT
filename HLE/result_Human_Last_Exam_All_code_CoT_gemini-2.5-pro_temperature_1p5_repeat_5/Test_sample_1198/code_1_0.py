import decimal

def demonstrate_construction(mod, steps, initial_M):
    """
    Demonstrates the constructive proof for the existence of 'a' for a few steps.
    
    Args:
        mod (int): The modulus (e.g., 2 or 3).
        steps (int): The number of steps (n) to demonstrate.
        initial_M (int): A starting integer value for floor(a).
    """
    # Set high precision for decimal calculations
    decimal.getcontext().prec = 100

    def target_residue(n, m):
        return n % m

    # --- Step n=1 ---
    # We need floor(a) ≡ 1 (mod m).
    # We choose a specific integer M for floor(a) and start with the interval [M, M+1).
    s1 = target_residue(1, mod)
    if initial_M % mod != s1:
        print(f"Warning: Initial M={initial_M} does not satisfy M % {mod} == {s1}.")
        # Adjust M to the next integer that satisfies the condition
        initial_M = initial_M - (initial_M % mod) + s1
        if initial_M % mod != s1:
            initial_M += mod
        print(f"Adjusted M to {initial_M} to satisfy the condition for n=1.")

    print(f"Starting with the condition floor(a) = {initial_M}, which means a is in [{initial_M}, {initial_M+1}).")
    # Use Decimal objects for high precision
    intervals = [[decimal.Decimal(initial_M), decimal.Decimal(initial_M + 1)]]
    print(f"n=1: Intervals for a: [[{intervals[0][0]}, {intervals[0][1]}]]\n")

    # --- Steps n > 1 ---
    for n in range(2, steps + 1):
        s_n = target_residue(n, mod)
        new_intervals = []
        print(f"--- Step n={n}, target residue for floor(a^{n}) is {s_n} ---")

        if not intervals:
            print("The set of possible 'a' became empty. This can happen if the initial M is not large enough.")
            break

        # Process each interval found in the previous step
        for L, R in intervals:
            # Calculate the range for a^n
            L_pow_n = L**n
            R_pow_n = R**n

            # Determine the range of integers to check for floor(a^n)
            k_start = int(L_pow_n.to_integral_value(rounding=decimal.ROUND_CEILING))
            k_end = int(R_pow_n.to_integral_value(rounding=decimal.ROUND_FLOOR))
            
            # Find all integers k in this range with the correct residue
            valid_k_found = []
            for k_int in range(k_start, k_end + 1):
                if k_int % mod == s_n:
                    k = decimal.Decimal(k_int)
                    # For a given k, a must be in [k^(1/n), (k+1)^(1/n))
                    new_L_cand = k.power(decimal.Decimal(1)/n)
                    new_R_cand = (k + 1).power(decimal.Decimal(1)/n)
                    
                    # The new interval for 'a' is the intersection of [L, R) and the interval for k
                    final_L = max(L, new_L_cand)
                    final_R = min(R, new_R_cand)

                    if final_L < final_R:
                        new_intervals.append([final_L, final_R])
                        valid_k_found.append(k_int)
            
            print(f"  For a in [{float(L):.12f}, {float(R):.12f}):")
            print(f"  Range for a^{n} is [{float(L_pow_n):.6f}, {float(R_pow_n):.6f}).")
            if valid_k_found:
                print(f"  Found valid integer(s) k: {valid_k_found}")
            else:
                 print(f"  No valid integer k with residue {s_n} found in this range.")


        intervals = new_intervals
        print(f"Result for n={n}: Found {len(intervals)} non-empty interval(s) for a.")
        # Print first few resulting intervals for brevity
        for i, (l,r) in enumerate(intervals[:5]):
             print(f"  Interval {i+1}: [{float(l):.12f}, {float(r):.12f})")
        if len(intervals) > 5:
            print(f"  ... and {len(intervals)-5} more intervals.")
        print("")

    if intervals:
        print("The process successfully found non-empty intervals at each step,")
        print("which supports the theoretical result that such a number 'a' exists.")


if __name__ == '__main__':
    print("This script demonstrates the existence of a real number 'a' for two related problems.")
    
    print("\n" + "="*60 + "\n")
    print("Part 1: Is there an a>0 such that floor(a^n) ≡ n mod 2 for all n>0?")
    print("="*60)
    demonstrate_construction(mod=2, steps=6, initial_M=3)
    
    print("\n" + "="*60 + "\n")
    print("Part 2: Is there an a>0 such that floor(a^n) ≡ n mod 3 for all n>0?")
    print("="*60)
    # A larger M and fewer steps are used, as the numbers grow faster.
    demonstrate_construction(mod=3, steps=5, initial_M=4)