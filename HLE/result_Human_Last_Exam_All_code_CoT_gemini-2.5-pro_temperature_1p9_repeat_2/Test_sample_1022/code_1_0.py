def compute_cohomology_dimension():
    """
    Computes the dimension of the degree 4 real cohomology group of T x T.

    This dimension corresponds to the dimension of the image of the comparison map
    c: H_b^4(T x T; R) -> H^4(T x T; R), providing a finite-dimensional
    answer to the posed problem.
    """
    
    # h_n stores the dimension of H^n(T; R), the nth real cohomology group of Thompson's group T.
    # Based on the known structure H*(T; Z) = Z[x_2, x_3]/(3x_3), with real coefficients the torsion part vanishes.
    h = {
        0: 1,  # dim H^0(T; R)
        1: 0,  # dim H^1(T; R)
        2: 1,  # dim H^2(T; R) = span(x_2)
        3: 0,  # dim H^3(T; R) is 0 because H^3(T;Z) is torsion
        4: 1,  # dim H^4(T; R) = span(x_2^2)
    }

    n = 4
    total_dim = 0
    
    print(f"Computing the dimension of H^{n}(T x T; R) using the Künneth formula:")
    print(f"dim H^{n}(T x T; R) = sum_{{p+q={n}}} dim(H^p(T; R)) * dim(H^q(T; R))")
    print("-" * 30)

    equation_parts = []
    for p in range(n + 1):
        q = n - p
        dim_p = h.get(p, 0)
        dim_q = h.get(q, 0)
        term_dim = dim_p * dim_q
        
        print(f"p={p}, q={q}: dim(H^{p}) * dim(H^{q}) = {dim_p} * {dim_q} = {term_dim}")
        
        equation_parts.append(str(term_dim))
        total_dim += term_dim
        
    final_equation = " + ".join(equation_parts)
    print("-" * 30)
    print(f"Total Dimension = {final_equation} = {total_dim}")

if __name__ == '__main__':
    compute_cohomology_dimension()
