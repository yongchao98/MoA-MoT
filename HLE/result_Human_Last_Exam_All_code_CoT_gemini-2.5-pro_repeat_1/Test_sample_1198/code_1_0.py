import math
from decimal import Decimal, getcontext

# Use Decimal for higher precision to handle small intervals
getcontext().prec = 50

def demonstrate_existence(mod, n_iterations, initial_k_start):
    """
    Demonstrates the existence of a real number 'a' such that
    floor(a^n) = n (mod m) for all n > 0 by showing that the
    interval of possible values for 'a' remains non-empty.

    Args:
        mod (int): The modulus (e.g., 2 or 3).
        n_iterations (int): How many steps (n) to run the demonstration for.
        initial_k_start (int): A starting integer for `a`, must be > mod.
    """
    print(f"--- Demonstrating existence for modulo {mod} ---")

    # Step n=1: Find the first integer k1 for floor(a)
    # We need floor(a) = k1, where k1 = 1 (mod m).
    # We choose a k1 > mod to ensure the range of a^n is large enough in subsequent steps.
    k1 = -1
    # Find the first integer >= initial_k_start that is 1 mod m
    for i in range(initial_k_start, initial_k_start + mod):
        if i % mod == 1:
            k1 = i
            break

    # Set the initial interval for 'a'
    c = Decimal(k1)
    d = Decimal(k1 + 1)
    print(f"n=1: We choose k1 = {k1} such that k1 > {mod} and k1 = 1 (mod {mod}).")
    print(f"     Initial interval for 'a' is [{c}, {d})")

    # Iterate for n > 1
    for n in range(2, n_iterations + 1):
        # Determine the required remainder for floor(a^n)
        required_rem = n % mod
        
        # Calculate the range for a^n given the current interval for 'a'
        cn = c**n
        dn = d**n
        
        # Search for a suitable integer k_n in [c^n, d^n)
        # We start our search from the first integer inside the range
        start_k_search = math.ceil(cn)
        found_k = False
        kn = -1
        
        # Iterate through candidate integers until one is found or the range is exhausted
        for k_candidate in range(start_k_search, math.ceil(dn)):
            if k_candidate % mod == required_rem:
                kn = k_candidate
                found_k = True
                break
        
        if not found_k:
            print(f"\nCould not find a suitable integer k{n} in range [{cn:.6f}, {dn:.6f}).")
            print("The constructive proof relies on this range being large enough.")
            return

        # A suitable k_n was found. Now refine the interval for 'a'.
        # The new condition is k_n <= a^n < k_n + 1, which means
        # k_n^(1/n) <= a < (k_n+1)^(1/n)
        new_c = Decimal(kn)**(Decimal(1)/Decimal(n))
        new_d = Decimal(kn + 1)**(Decimal(1)/Decimal(n))

        # The updated interval for 'a' is the intersection of the old and new intervals
        c = max(c, new_c)
        d = min(d, new_d)

        if c >= d:
            print(f"\nError: Interval became empty at n={n}. This shouldn't happen with high enough precision.")
            return
            
        print(f"n={n}: Required remainder is {required_rem}. Found k{n} = {kn}.")
        print(f"     New interval for 'a' is [{c:.15f}, {d:.15f})")
        print(f"     Interval width: {d-c:.2E}")

# --- Main execution ---
# Case: modulo 2
demonstrate_existence(mod=2, n_iterations=10, initial_k_start=3)

print("\n" + "="*60 + "\n")

# Case: modulo 3
demonstrate_existence(mod=3, n_iterations=10, initial_k_start=4)