import math

def find_a_intervals(m, N, initial_intervals):
    """
    Tries to find intervals for a real number 'a' satisfying the condition
    floor(a^n) % m == n % m for n from 1 to N.

    Args:
        m (int): The modulus (2 or 3).
        N (int): The number of steps to iterate.
        initial_intervals (list of tuples): The starting range for 'a'.
    """
    print(f"--- Searching for a number 'a' with floor(a^n) % {m} == n % {m} ---")
    current_intervals = initial_intervals
    
    for n in range(1, N + 1):
        next_intervals = []
        target_rem = n % m
        
        if not current_intervals:
            print(f"n={n}: Search failed. No possible intervals for 'a' remain.")
            return []

        # Iterate through all current valid intervals for 'a'
        for L, R in current_intervals:
            # For a given interval (L, R), the value of k = floor(a^n)
            # can range from floor(L^n) to ceil(R^n).
            # We add a buffer to k_end for safety with floating point math.
            k_start = math.floor(L**n)
            k_end = math.ceil(R**n) + 1
            
            # We are looking for an integer k in this range with k % m == target_rem
            for k in range(k_start, k_end):
                if k > 0 and k % m == target_rem:
                    # The condition floor(a^n) = k implies k <= a^n < k+1,
                    # which means a is in the interval [k^(1/n), (k+1)^(1/n)).
                    nL = k**(1/n)
                    nR = (k+1)**(1/n)
                    
                    # We intersect this new interval with our current interval [L, R)
                    # to find the updated valid range for 'a'.
                    intersect_L = max(L, nL)
                    intersect_R = min(R, nR)
                    
                    # If the intersection is a valid interval, add it to our list
                    # for the next iteration.
                    if intersect_L < intersect_R:
                        next_intervals.append((intersect_L, intersect_R))
        
        current_intervals = next_intervals
        
        print(f"Step n={n}:")
        if current_intervals:
            print(f"  Condition: floor(a^{n}) % {m} == {target_rem}")
            # The number of intervals can grow large, so we just show a summary.
            print(f"  Found {len(current_intervals)} possible interval(s) for 'a'.")
            # To keep it simple, we print the first one found.
            first_interval = current_intervals[0]
            print(f"  An example interval is: [{first_interval[0]:.9f}, {first_interval[1]:.9f})")
        else:
            print(f"  No interval found. It's not possible to satisfy the condition for n={n}.")
            break
            
    return current_intervals

if __name__ == '__main__':
    # --- Case 1: Modulo 2 ---
    # For n=1, we need floor(a) % 2 == 1. Let's start with a in [3.0, 4.0).
    initial_intervals_m2 = [(3.0, 4.0)]
    find_a_intervals(m=2, N=10, initial_intervals=initial_intervals_m2)

    print("\n" + "="*50 + "\n")

    # --- Case 2: Modulo 3 ---
    # For n=1, we need floor(a) % 3 == 1. Let's start with a in [4.0, 5.0).
    initial_intervals_m3 = [(4.0, 5.0)]
    find_a_intervals(m=3, N=10, initial_intervals=initial_intervals_m3)