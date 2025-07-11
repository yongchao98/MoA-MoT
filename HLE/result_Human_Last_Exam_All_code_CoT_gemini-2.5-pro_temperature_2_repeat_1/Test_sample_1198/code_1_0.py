import math

def get_new_intervals(current_intervals, n, q):
    """
    Given a list of intervals for 'a', this function calculates the new,
    more constrained list of intervals based on the condition for step 'n'.
    """
    next_intervals = []
    print(f"--- Step n={n}, q={q} ---")
    
    if not current_intervals:
        print("No possible intervals for 'a' from previous step.")
        return []

    print(f"Current interval(s) for a: {current_intervals}")

    # For each interval in our list of possibilities for 'a'
    for l, r in current_intervals:
        # Determine the range for a^n
        l_pow_n = l**n
        r_pow_n = r**n
        print(f"  For a in [{l:.4f}, {r:.4f}), a^{n} is in [{l_pow_n:.4f}, {r_pow_n:.4f})")

        # Determine the possible integer values for floor(a^n)
        min_floor = math.ceil(l_pow_n)
        max_floor = math.floor(r_pow_n - 1e-9) # A small epsilon to handle right boundary
        
        possible_floors = range(min_floor, max_floor + 1)
        if not possible_floors:
            print(f"  Range for a^{n} contains no integers. Cannot continue.")
            continue
            
        print(f"  Possible floors for a^{n}: {list(possible_floors)}")

        # Find which of these floors satisfy the congruence
        target_residue = n % q
        valid_floors = [m for m in possible_floors if m % q == target_residue]
        print(f"  We need floor(a^{n}) % {q} == {target_residue}.")
        print(f"  Valid floors: {valid_floors}")
        
        if not valid_floors:
            print(f"  No valid floor found. This path fails.")
            continue

        # For each valid floor, find the corresponding interval for 'a' and intersect it with the current one.
        for m in valid_floors:
            # New constraint on 'a' is m <= a^n < m+1, which is m^(1/n) <= a < (m+1)^(1/n)
            new_l = m**(1/n)
            new_r = (m+1)**(1/n)
            
            # Intersect [new_l, new_r) with the current interval [l, r)
            intersect_l = max(l, new_l)
            intersect_r = min(r, new_r)

            if intersect_l < intersect_r:
                next_intervals.append((intersect_l, intersect_r))

    print(f"New interval(s) for a: {next_intervals}")
    return next_intervals

def main():
    """
    Main function to run the demonstration for modulo 3.
    """
    q = 3
    # A list of tuples, representing the possible disjoint intervals for 'a'.
    # Start with the condition for n=1: floor(a) = 1 (mod 3).
    # We choose the simplest case: floor(a) = 1, so 1 <= a < 2.
    # Note: We could have chosen floor(a) = 4, etc. which would lead to other branches.
    intervals_a = [(1.0, 2.0)]
    
    # We could have also chosen floor(a)=4 --> a in [4,5).
    # Try running the code with intervals_a = [(4.0, 5.0)] to see another path.

    for n in range(1, 7):
        intervals_a = get_new_intervals(intervals_a, n, q)
        if not intervals_a:
            print("\nResult: The set of possible values for 'a' became empty.")
            print(f"This demonstrates that it is impossible to satisfy the condition for all n.")
            break
        print("\n")

if __name__ == '__main__':
    main()
