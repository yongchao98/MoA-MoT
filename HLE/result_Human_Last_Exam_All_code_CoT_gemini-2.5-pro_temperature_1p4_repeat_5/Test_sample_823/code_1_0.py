def min_vertices_for_induced_matching(k, d):
    """
    Calculates the minimum number of vertices n required to guarantee an
    induced matching of size at least k in a graph with maximum degree d.

    The argument is based on the inequality: size_of_matching >= n / (2*d + 2)
    To guarantee a matching of size k, we need n / (2*d + 2) > k - 1,
    or more simply, ensure the lower bound is at least k. So we need n >= k * (2*d + 2).
    """
    # To guarantee a matching of size k, we need n to be large enough.
    # n_min is the smallest integer n such that floor(n / (2d+2)) >= k
    # which is satisfied by n_min = k * (2d + 2).
    # Any graph with more vertices will also satisfy this.
    n_min = k * (2 * d + 2)
    guaranteed_size = n_min / (2 * d + 2)

    print(f"Given a maximum degree d = {d}.")
    print(f"To guarantee an induced matching of size at least k = {k}:")
    print(f"A graph must have at least n = {n_min} vertices.")
    print(f"Because the size of the induced matching is guaranteed to be at least:")
    print(f"n / (2*d + 2) = {n_min} / (2*{d} + 2) = {int(guaranteed_size)}")


# Example from the problem description
# Let's say d is some constant, e.g., d=4
# And we want to find a graph with an induced matching of size k=100.
d_constant = 4
k_target = 100
min_vertices_for_induced_matching(k_target, d_constant)