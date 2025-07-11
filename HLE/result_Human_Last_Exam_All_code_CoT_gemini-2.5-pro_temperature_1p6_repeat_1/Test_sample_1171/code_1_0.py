def main():
    """
    Finds the values of k in {1, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.
    """

    # k values for which pi_k(S^4) @ Q is non-zero.
    # For S^{2n}, non-zero only at k=2n and k=4n-1. Here n=2.
    s4_nonzero_k = {4, 2 * 4 - 1}

    # k values for which pi_k(CP^2) @ Q is non-zero.
    # For CP^n, non-zero only at k=2 and k=2n+1. Here n=2.
    cp2_nonzero_k = {2, 2 * 2 + 1}

    vanishing_k = []
    for k in range(1, 10):
        # The rational homotopy group of the wedge sum vanishes iff
        # the groups for both components vanish.
        is_s4_zero = k not in s4_nonzero_k
        is_cp2_zero = k not in cp2_nonzero_k

        if is_s4_zero and is_cp2_zero:
            vanishing_k.append(k)

    # Print the final list of k values as a comma-separated string.
    print(*vanishing_k, sep=',')

if __name__ == "__main__":
    main()