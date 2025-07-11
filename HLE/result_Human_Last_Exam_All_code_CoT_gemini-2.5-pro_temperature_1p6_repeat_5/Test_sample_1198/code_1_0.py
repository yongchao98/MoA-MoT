import math
import decimal

def find_a_mod_2(n_max=20, initial_k1=3):
    """
    This function demonstrates the constructive proof for the existence of a real number `a`
    such that floor(a^n) = n (mod 2) for all n > 0.
    It constructs a series of nested intervals for `a`.
    """
    # Use a high precision context for calculations to handle the small intervals
    decimal.getcontext().prec = 100

    # Step 1: n=1. We need floor(a) to be odd.
    # We choose an initial odd integer, e.g., initial_k1.
    # This defines the initial interval for `a` as [initial_k1, initial_k1 + 1).
    k1 = initial_k1
    lower_bound = decimal.Decimal(k1)
    upper_bound = decimal.Decimal(k1 + 1)
    
    print(f"Plan: Construct nested intervals for a, ensuring floor(a^n) has the right parity.\n")
    print(f"Step n=1: Required parity is 1. We choose k_1 = floor(a) = {k1}.")
    print(f"Resulting interval for a: [{lower_bound}, {upper_bound})\n")

    # Loop for n = 2, 3, ...
    for n in range(2, n_max + 1):
        required_parity = n % 2
        
        # Calculate the range for a^n based on the current interval for a.
        an_lower_bound = lower_bound ** n
        an_upper_bound = upper_bound ** n
        
        # We need to find an integer k_n in [an_lower_bound, an_upper_bound)
        # such that k_n has the required parity.
        
        # Start searching for a suitable integer k_n from the lower end of the range.
        # We search through integers m in the range [floor(an_lower_bound), ceil(an_upper_bound)]
        start_m = math.floor(an_lower_bound)
        end_m = math.ceil(an_upper_bound)
        
        k_n = None
        for m_candidate in range(start_m, end_m + 1):
            if m_candidate % 2 == required_parity:
                # This integer has the correct parity. Now, will it produce a valid, non-empty interval for a?
                # The interval for a would be [m_candidate^(1/n), (m_candidate+1)^(1/n)).
                # We need to check if its intersection with the current interval [lower_bound, upper_bound) is non-empty.
                
                # Calculate the potential new interval bounds for a.
                potential_lower = decimal.Decimal(m_candidate) ** (decimal.Decimal(1)/n)
                potential_upper = decimal.Decimal(m_candidate + 1) ** (decimal.Decimal(1)/n)
                
                # Check for non-empty intersection.
                if max(lower_bound, potential_lower) < min(upper_bound, potential_upper):
                    k_n = m_candidate
                    break # We found a suitable k_n, so we'll use it.
        
        print(f"Step n={n}: Required parity is {n % 2}.")
        if k_n is None:
            print("Failed to find a suitable integer k_n. The interval has become empty.")
            # This case shouldn't happen based on the theory, but it's a safeguard.
            return

        print(f"Range for a^{n} is approx. [{float(an_lower_bound):.4f}, {float(an_upper_bound):.4f}). We choose k_{n} = {k_n}.")

        # Update the interval for `a` by intersecting it with the new constraint.
        new_lower = decimal.Decimal(k_n) ** (decimal.Decimal(1)/n)
        new_upper = decimal.Decimal(k_n + 1) ** (decimal.Decimal(1)/n)
        lower_bound = max(lower_bound, new_lower)
        upper_bound = min(upper_bound, new_upper)
        
        print(f"Resulting interval for a: [{float(lower_bound):.15f}, {float(upper_bound):.15f})\n")

find_a_mod_2()