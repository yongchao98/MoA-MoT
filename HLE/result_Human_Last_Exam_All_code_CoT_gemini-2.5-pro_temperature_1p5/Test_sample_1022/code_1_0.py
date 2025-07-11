def solve_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.
    
    This function uses the Künneth formula for bounded cohomology and the fact
    that the bounded cohomology of Thompson's group T is trivial in positive degrees.
    """
    
    # The degree of the cohomology group we want to compute.
    n = 4
    
    # Known dimensions of the bounded cohomology groups of Thompson's group T.
    # dim H_b^0(T, R) = 1
    # dim H_b^k(T, R) = 0 for k > 0
    # We need dimensions for degrees 0 through n.
    dim_Hb_T = [0] * (n + 1)
    if n >= 0:
        dim_Hb_T[0] = 1
        
    print(f"This script computes the dimension of the degree {n} bounded cohomology group of T x T with trivial real coefficients.")
    print("The calculation is based on the Künneth formula for bounded cohomology:")
    print(f"dim H_b^{n}(T x T) = sum_{{p+q={n}}} (dim H_b^p(T) * dim H_b^q(T))")
    print("\nKnown dimensions for Thompson's group T:")
    for i in range(n + 1):
        print(f"dim H_b^{i}(T) = {dim_Hb_T[i]}")

    total_dimension = 0
    equation_parts = []
    
    # Loop through all possible values of p from 0 to n.
    for p in range(n + 1):
        q = n - p
        
        # Get dimensions for H_b^p(T) and H_b^q(T)
        dim_p = dim_Hb_T[p]
        dim_q = dim_Hb_T[q]
        
        # Calculate the dimension of the tensor product term
        term_dimension = dim_p * dim_q
        total_dimension += term_dimension
        
        # Store the string representation for the final equation
        equation_parts.append(f"({dim_p} * {dim_q})")
        
    # Construct and print the final equation with numerical values.
    equation_str = " + ".join(equation_parts)
    print("\nApplying the formula for n=4:")
    print(f"dim H_b^4(T x T) = {equation_str}")
    
    # Print the sum of the calculated terms
    term_values = [str(eval(part)) for part in equation_parts]
    term_sum_str = " + ".join(term_values)
    print(f"                 = {term_sum_str}")
    
    # Print the final result
    print(f"                 = {total_dimension}")
    print("\n-------------------------------------------------------------")
    print(f"The dimension of the degree {n} bounded cohomology group of T x T is {total_dimension}.")
    print("-------------------------------------------------------------")

if __name__ == "__main__":
    solve_cohomology_dimension()
