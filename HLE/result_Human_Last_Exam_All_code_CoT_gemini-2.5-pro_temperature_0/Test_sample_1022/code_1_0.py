def dim_Hb_T(k):
    """
    Computes the dimension of the k-th bounded cohomology group of Thompson's group T
    with trivial real coefficients, dim H_b^k(T, R).

    According to a theorem by Burger and Monod, the bounded cohomology ring
    H_b^*(T, R) is a polynomial algebra R[Ω] on a single generator Ω of degree 3.
    """
    if not isinstance(k, int) or k < 0:
        return 0
    # dim H_b^0(G, R) is always 1 for any group G.
    if k == 0:
        return 1
    # For k > 0, the dimension is 1 if k is a multiple of 3, and 0 otherwise.
    if k > 0 and k % 3 == 0:
        return 1
    return 0

def solve_and_print():
    """
    Computes and prints the dimension of the degree 4 bounded cohomology group
    of the direct product T x T.
    """
    n = 4
    print(f"We aim to compute the dimension of H_b^{n}(T x T, R).")
    print("Based on the Künneth formula for bounded cohomology, the dimension is given by the sum:")
    print(f"dim H_b^{n}(T x T, R) = sum_{{i+j={n}}} dim H_b^i(T, R) * dim H_b^j(T, R)")
    print("\nFirst, let's list the necessary dimensions for H_b^k(T, R):")
    for k in range(n + 1):
        print(f"  dim H_b^{k}(T, R) = {dim_Hb_T(k)}")

    print(f"\nNow, we compute each term in the sum for n={n}:")
    
    terms = []
    for i in range(n + 1):
        j = n - i
        dim_i = dim_Hb_T(i)
        dim_j = dim_Hb_T(j)
        term_product = dim_i * dim_j
        terms.append(term_product)
        print(f"  Term (i={i}, j={j}): dim H_b^{i}(T, R) * dim H_b^{j}(T, R) = {dim_i} * {dim_j} = {term_product}")

    total_dimension = sum(terms)
    
    # Format the final equation as requested
    equation_str = " + ".join(map(str, terms))
    
    print("\nThe final calculation is the sum of these terms:")
    print(f"{equation_str} = {total_dimension}")
    
    print(f"\nThus, the dimension of the degree {n} bounded cohomology group of T x T is {total_dimension}.")

if __name__ == "__main__":
    solve_and_print()