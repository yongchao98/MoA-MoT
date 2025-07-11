import math

def solve_for_modulus(modulus):
    """
    Demonstrates the search for a number 'a' for a given modulus.
    It follows a constructive proof method by narrowing down intervals for 'a'.
    """
    print(f"--- Testing for modulus {modulus} ---")
    
    # For n=1, we need floor(a) % modulus == 1.
    # We choose the smallest possible integer value for floor(a), which is 1.
    # This means we can start with the assumption that 'a' is in the interval [1.0, 2.0).
    intervals = [(1.0, 2.0)]
    print(f"Step n=1: We need floor(a) % {modulus} == 1.")
    print(f"Assuming floor(a)=1, the interval for 'a' is [1.0, 2.0).")

    n_max = 7 # Check up to n=7, which is sufficient to illustrate the process.

    for n in range(2, n_max + 1):
        # The remainder we need for floor(a^n)
        required_remainder = n % modulus
        
        print(f"\nStep n={n}: We need floor(a^{n}) % {modulus} == {required_remainder}.")

        next_intervals = []
        for L, R in intervals:
            # For an interval [L, R) for 'a', a^n is in [L^n, R^n).
            lower_bound_an = L**n
            upper_bound_an = R**n

            # We need to find integers 'm' in [lower_bound_an, upper_bound_an)
            # such that m % modulus == required_remainder.
            
            # Start checking from the first integer >= lower_bound_an
            m_start = math.ceil(lower_bound_an)

            # End checking at the first integer < upper_bound_an
            m_end = math.floor(upper_bound_an)
            if m_end == upper_bound_an:
                m_end -= 1

            for m in range(m_start, m_end + 1):
                if m % modulus == required_remainder:
                    # If we find such an m, it's a candidate for floor(a^n).
                    # This implies m <= a^n < m+1.
                    # Taking the n-th root gives an interval for 'a'.
                    new_L = m**(1/n)
                    new_R = (m+1)**(1/n)

                    # The new interval for 'a' must also be within [L, R).
                    # So we take the intersection.
                    intersect_L = max(L, new_L)
                    intersect_R = min(R, new_R)

                    # If the intersection is valid (not empty), add it to our list.
                    if intersect_L < intersect_R:
                        next_intervals.append((intersect_L, intersect_R))

        intervals = next_intervals
        
        if not intervals:
            print("The set of possible intervals became empty.")
            print(f"This suggests no such 'a' exists (within the initial range).")
            return

        print(f"The set of possible 'a' values is now the union of {len(intervals)} interval(s):")
        # Print a few examples for brevity
        for i in range(min(len(intervals), 3)):
             print(f"  ... e.g., [{intervals[i][0]:.8f}, {intervals[i][1]:.8f})")

    print("\nAfter several steps, the set of intervals is still non-empty.")
    print("This provides numerical evidence that such a number 'a' likely exists.")

# Run the demonstration for both cases
solve_for_modulus(2)
print("\n" + "="*40 + "\n")
solve_for_modulus(3)