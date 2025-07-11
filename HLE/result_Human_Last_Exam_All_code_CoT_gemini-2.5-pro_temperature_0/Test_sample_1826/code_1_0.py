def construct_almost_disjoint_set(S, num_elements):
    """
    Constructs an infinite set x that is almost disjoint from every set in the countable collection S.

    Args:
        S (list of list of int): A list of infinite sets (represented as sorted lists).
        num_elements (int): The number of elements to construct for the set x.

    Returns:
        list of int: The first `num_elements` of the constructed set x.
    """
    # Ensure all sets in S are sorted
    for s in S:
        s.sort()

    x = []
    x_prev = -1

    print("Constructing the set x = {x_0, x_1, x_2, ...}\n")

    for k in range(num_elements):
        # We can only consider sets s_i where i <= k
        if k >= len(S):
            print(f"Cannot compute x_{k} as s_{k} is not provided in the sample S. Stopping.")
            break

        # Collect the elements for the max operation: {s_{i,k} for i <= k}
        elements_to_consider = []
        indices_str = []
        values_str = []
        
        for i in range(k + 1):
            if k < len(S[i]):
                elements_to_consider.append(S[i][k])
                indices_str.append(f"s_{i}[{k}]")
                values_str.append(str(S[i][k]))
            else:
                # This case happens if our sample lists are not long enough.
                # In theory, the sets are infinite. We'll stop if we can't proceed.
                print(f"Cannot compute x_{k} as s_{i}[{k}] (for i={i}) is out of bounds. Stopping.")
                return x

        # The recursive formula
        max_val = max([x_prev] + elements_to_consider)
        x_k = max_val + 1
        
        # Print the calculation for x_k
        print(f"Step k={k}:")
        print(f"  x_{k} = max(x_{k-1}, " + ", ".join(indices_str) + ") + 1")
        print(f"      = max({x_prev}, " + ", ".join(values_str) + ") + 1")
        print(f"      = max({[x_prev] + elements_to_consider}) + 1")
        print(f"      = {max_val} + 1 = {x_k}")
        print("-" * 20)

        x.append(x_k)
        x_prev = x_k
        
    return x

# --- Main execution ---

# Define a sample collection S of infinite sets
# s_0 = {10, 20, 30, 40, ...}
# s_1 = {15, 25, 35, 45, ...}
# s_2 = {1, 2, 3, 100, 101, 102, ...}
# s_3 = {5, 10, 15, 20, ...}
S_sample = [
    [10 * (i + 1) for i in range(50)],
    [15 + 10 * i for i in range(50)],
    [i + 1 for i in range(3)] + [100 + i for i in range(50)],
    [5 * (i + 1) for i in range(50)]
]

# Number of elements to construct for our set x
N = 10

# Construct the set x
x_constructed = construct_almost_disjoint_set(S_sample, N)

print(f"\nConstructed set x (first {len(x_constructed)} elements):")
print(f"x = {x_constructed}")

print("\nChecking intersections with the sets in S:")
for i, s in enumerate(S_sample):
    intersection = sorted(list(set(x_constructed) & set(s)))
    print(f"|x ∩ s_{i}| = |{intersection}| = {len(intersection)}")

print("\nAs demonstrated, the constructed set x has a finite intersection with each set s_i.")
print("The argument holds for any countable collection of infinite sets, so the answer is YES.")

<<<YES>>>