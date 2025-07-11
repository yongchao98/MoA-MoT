def find_almost_disjoint_set(s_list, num_elements):
    """
    Given a countable family of infinite sets S (represented by a list of lists),
    this function constructs the first `num_elements` of an infinite set x
    that is almost disjoint from every set in S.

    The method is based on constructing a "dominating function" g,
    and the set x is the range of g.
    """
    # s_list represents the family S = {s_0, s_1, s_2, ...}
    # The lists act as the enumerating functions for the sets.
    # f_n(k) is equivalent to s_list[n][k].

    g = []
    print("Construction of the dominating function g:")
    print("-" * 40)
    print(f"{'k':>3} | {'h(k)':>8} | {'g(k)':>8} | Computation")
    print("-" * 40)

    for k in range(num_elements):
        # Step 1: Calculate h(k) = max(f_0(k), ..., f_k(k)) + 1
        # This requires the first k+1 lists to have at least k+1 elements.
        # We handle potential index errors for this finite demonstration.
        try:
            max_val = -1
            f_k_values = []
            # We can only go up to min(k, len(s_list) - 1) for the set index
            limit = min(k, len(s_list) - 1)
            for n in range(limit + 1):
                val = s_list[n][k]
                f_k_values.append(f"f_{n}({k})={val}")
                if val > max_val:
                    max_val = val
            h_k = max_val + 1
            h_k_comp = f"max({', '.join(f_k_values)})+1 = {h_k}"
        except IndexError:
            print(f"Stopping at k={k} due to insufficient length of input lists.")
            break

        # Step 2: Calculate g(k) ensuring it is strictly increasing
        if k == 0:
            g_k = h_k
            g_k_comp = f"g(0) = h(0) = {g_k}"
        else:
            prev_g = g[k-1]
            g_k = max(h_k, prev_g + 1)
            g_k_comp = f"g({k}) = max(h({k}), g({k-1})+1) = max({h_k}, {prev_g}+1) = {g_k}"

        g.append(g_k)
        print(f"{k:>3} | {h_k:>8} | {g_k:>8} | {g_k_comp}")

    print("-" * 40)
    return g

# --- Example Usage ---
# Let's define a sample family of sets S = {s_0, s_1, s_2}
# s_0: Prime numbers
s_0 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
# s_1: Perfect squares
s_1 = [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196]
# s_2: Fibonacci numbers (unique and sorted)
s_2 = [0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610]

S = [s_0, s_1, s_2]
num_x_elements = 10
x_set = find_almost_disjoint_set(S, num_x_elements)

print("\n--- Results ---")
print("The family of sets S (first 15 elements):")
for i, s in enumerate(S):
    print(f"s_{i}: {s}")

print(f"\nThe constructed set x (first {num_x_elements} elements):")
print(f"x = {x_set}")

print("\nAnalysis of intersections with the original sets:")
for i, s in enumerate(S):
    # Use the full provided list for s for a better check
    intersection = sorted(list(set(x_set) & set(s)))
    print(f"x ∩ s_{i} = {intersection} (size: {len(intersection)})")