def vogel_braid_index_bound():
    """
    Calculates an upper bound for the braid index of the three-twist knot (5_2)
    using Vogel's algorithm.
    """
    # The Alexander polynomial for the three-twist knot is 2t - 3 + 2t^-1.
    # We list the coefficients in order of increasing powers (t^-1, t^0, t^1).
    coeffs = [2, -3, 2]
    min_degree = -1

    print("1. The Alexander polynomial for the three-twist knot is Δ(t) = 2t^{-1} - 3t^{0} + 2t^{1}.")
    print(f"The coefficients are: {coeffs}")
    print("\n2. We compute the partial sums (S_k) of the coefficients.")

    # S_k = 0 for k < min_degree
    partial_sums = [0]
    running_sum = 0
    print("   S_k = 0 for k < -1")
    
    # Calculate sums for non-zero coefficients
    for i, c in enumerate(coeffs):
        running_sum += c
        partial_sums.append(running_sum)
        print(f"   S_{min_degree + i} = {running_sum}")
        
    print(f"\nThe set of distinct partial sums is: {sorted(list(set(partial_sums)))}")

    # Find the maximum and minimum of the partial sums
    max_s = max(partial_sums)
    min_s = min(partial_sums)

    print("\n3. Find the maximum and minimum of the partial sums:")
    print(f"   Maximum value = {max_s}")
    print(f"   Minimum value = {min_s}")

    # Calculate the upper bound using Vogel's formula
    upper_bound = 1 + max_s - min_s

    print("\n4. The upper bound for the braid index is given by: 1 + max(S_k) - min(S_k)")
    print(f"\nFinal calculation: {1} + {max_s} - ({min_s}) = {upper_bound}")


if __name__ == "__main__":
    vogel_braid_index_bound()