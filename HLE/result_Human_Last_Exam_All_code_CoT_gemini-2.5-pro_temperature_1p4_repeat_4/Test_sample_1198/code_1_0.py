import math

def investigate_existence(q, n_max, initial_interval):
    """
    This function investigates if a real number 'a' can be found such that
    floor(a^n) = n (mod q) for all n > 0.
    It demonstrates the search process by iteratively narrowing down possible
    intervals for 'a'. For q=3, this process fails on the tested path.
    """
    
    # intervals are stored as a list of tuples (low, high), representing [low, high).
    intervals = [initial_interval]
    
    print(f"--- Investigation for modulus q = {q} ---")
    print(f"Starting search with 'a' in the interval {initial_interval}\n")

    for n in range(1, n_max + 1):
        rem = n % q
        
        print(f"Step n={n}:")
        print(f"  Condition: floor(a^{n}) % {q} == {rem}")
        
        if not intervals:
            print("  The list of valid intervals is empty. Search cannot continue on this path.")
            break
        
        print(f"  Current interval(s) for 'a': {[(f'{l:.8f}', f'{h:.8f}') for l, h in intervals]}")

        next_intervals = []
        
        for low, high in intervals:
            # Range of a^n
            a_n_low = low**n
            a_n_high = high**n
            
            # Find possible integer values for floor(a^n)
            m_min = math.ceil(a_n_low)
            # Use a small epsilon to handle boundary cases where a_n_high is an exact integer
            m_max = math.floor(a_n_high - 1e-9) 

            # Filter integers by the modulo condition
            valid_m = [m for m in range(m_min, m_max + 1) if m % q == rem]
            
            if not valid_m:
                print(f"  -> Path failed for interval [{low:.8f}, {high:.8f}):")
                print(f"     a^{n} is in [{a_n_low:.8f}, {a_n_high:.8f})")
                print(f"     The integers in this range are {list(range(m_min, m_max + 1))}, but none satisfy the condition.")
                continue

            # Calculate new intervals for 'a' and intersect with the current interval [low, high)
            for m in valid_m:
                new_low = m**(1/n)
                new_high = (m+1)**(1/n)
                
                intersect_low = max(low, new_low)
                intersect_high = min(high, new_high)
                
                if intersect_low < intersect_high:
                    next_intervals.append((intersect_low, intersect_high))
        
        intervals = next_intervals
        print("-" * 25)
        
        if not intervals:
            print(f"\nResult: After step n={n}, the set of possible values for 'a' has become empty.")
            print("This search path leads to a contradiction, demonstrating that no such 'a' exists on this path.")
            return

    if intervals:
      print(f"\nAfter {n_max} steps, the following interval(s) for 'a' remain valid on this path:")
      print(intervals)


# It is known that for q=3, no such 'a' exists. We demonstrate the failure of one search path.
investigate_existence(q=3, n_max=5, initial_interval=(1.0, 2.0))
