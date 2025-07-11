import sys

def solve_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group
    of T x T, where T is Thompson's group.
    """
    
    # We want to compute the dimension of H_b^4(T x T).
    n = 4
    
    # We use the known dimensions of the bounded cohomology groups of T.
    # dim H_b^k(T; R) for k=0, 1, 2, 3, 4
    # For k > 4, the dimension is also 0.
    dim_Hb_T = {
        0: 1,
        1: 0,
        2: 1,
        3: 1,
        4: 0,
    }
    
    print("This script computes the dimension of the degree 4 bounded cohomology group of T x T.")
    print("The calculation is based on the Künneth-type formula for bounded cohomology:")
    print(f"dim H_b^{n}(T x T) = sum_{{p+q={n}}} dim H_b^p(T) * dim H_b^q(T)\n")

    total_dim = 0
    equation_terms = []
    
    # Iterate through p from 0 to n. q is determined by p + q = n.
    for p in range(n + 1):
        q = n - p
        
        dim_p = dim_Hb_T.get(p, 0)
        dim_q = dim_Hb_T.get(q, 0)
        
        term = dim_p * dim_q
        total_dim += term
        
        # Store the term as a string for the final equation
        equation_terms.append(f"{dim_p}*{dim_q}")
        
    # Build the final equation string from the terms
    equation_string = " + ".join(equation_terms)
    
    print("The final computation is:")
    print(f"{equation_string} = {total_dim}")

if __name__ == "__main__":
    solve_cohomology_dimension()
