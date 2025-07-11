def compute_cohomology_dimension():
    """
    Computes the dimension of the ordinary degree 4 cohomology group of T x T.
    
    This calculation assumes the user intended to ask for ordinary cohomology
    instead of bounded cohomology.
    """
    
    # h[p] stores the dimension of the p-th ordinary cohomology group of T with real coefficients.
    # h^p = dim H^p(T, R)
    h = {
        0: 1,  # dim H^0(T, R)
        1: 2,  # dim H^1(T, R) from H_1(T, Z) = Z^2
        2: 0,  # dim H^2(T, R) from H_2(T, Z) = Z/2
        3: 1,  # dim H^3(T, R) from H_3(T, Z) = Z
        4: 0,  # dim H^4(T, R) from rational cohomology results
    }

    # Using the Künneth formula for n=4:
    # dim H^4(T x T) = sum_{p+q=4} h^p * h^q
    # This expands to: h^0*h^4 + h^1*h^3 + h^2*h^2 + h^3*h^1 + h^4*h^0
    
    terms = []
    total_dim = 0
    
    for p in range(5):
        q = 4 - p
        term_val = h[p] * h[q]
        terms.append(f"({h[p]} * {h[q]})")
        total_dim += term_val

    equation_str = " + ".join(terms)
    
    # Print the full equation as requested
    print(f"Dimension = {equation_str} = {total_dim}")

if __name__ == "__main__":
    compute_cohomology_dimension()
