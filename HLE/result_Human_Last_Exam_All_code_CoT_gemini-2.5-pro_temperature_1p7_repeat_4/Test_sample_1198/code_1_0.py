import math
import decimal

# Set precision for Decimal calculations for high accuracy
decimal.getcontext().prec = 50

def demonstrate_construction(mod, k1_start, steps=15):
    """
    Demonstrates the constructive proof for the existence of 'a' by building
    a sequence of nested intervals where 'a' must lie.
    """
    print(f"\n--- Demonstration for modulo {mod} ---")
    
    # We choose k1 such that k1 > mod and k1 % mod == 1 to guarantee success.
    if k1_start % mod != 1 or k1_start <= mod:
        print(f"Error: k1_start must be > {mod} and satisfy k1_start % {mod} == 1.")
        return

    # Step 1: Initialize with k1
    k_n = decimal.Decimal(k1_start)
    
    # The first interval for 'a' is [k1, k1+1)
    L, R = k_n, k_n + 1
    print(f"n = 1: Chose k1 = {k_n}. Interval for 'a': [{float(L):.12f}, {float(R):.12f})")

    # Step 2: Iterate for n = 2, 3, ... to refine the interval
    for n in range(2, steps + 1):
        # For the next step n, k_n must be an integer within the range [L^n, R^n).
        k_min_val = L**n
        k_max_val = R**n

        k_min = math.ceil(k_min_val)
        target_rem = n % mod
        
        # Find the first integer >= k_min that satisfies the congruence k_n % mod == n % mod
        first_candidate = k_min
        if first_candidate % mod != target_rem:
             first_candidate += (target_rem - (first_candidate % mod) + mod) % mod

        # Check if this candidate is within our allowed range
        if first_candidate < k_max_val:
            k_n = decimal.Decimal(first_candidate)
            # Update the interval [L, R) for 'a'. The new interval is [k_n^(1/n), (k_n+1)^(1/n))
            L = k_n**(decimal.Decimal(1)/n)
            R = (k_n + 1)**(decimal.Decimal(1)/n)
            print(f"n = {n:2}: Found k{n:<5} = {k_n}. New interval: [{float(L):.12f}, {float(R):.12f})")
        else:
            # This case shouldn't be reached with our choice of k1.
            print(f"n = {n}: Construction failed. No suitable k{n} found.")
            print(f"       Search range for k{n}: [{float(k_min_val):.4f}, {float(k_max_val):.4f})")
            print(f"       Interval length: {float(k_max_val - k_min_val):.4f}")
            return

# For modulo 2, the condition is floor(a^n) % 2 == n % 2.
# We need a > 2. We start with k1 > 2 and k1 % 2 == 1, so we choose k1=3.
demonstrate_construction(mod=2, k1_start=3)

# For modulo 3, the condition is floor(a^n) % 3 == n % 3.
# We need a > 3. We start with k1 > 3 and k1 % 3 == 1, so we choose k1=4.
demonstrate_construction(mod=3, k1_start=4)
