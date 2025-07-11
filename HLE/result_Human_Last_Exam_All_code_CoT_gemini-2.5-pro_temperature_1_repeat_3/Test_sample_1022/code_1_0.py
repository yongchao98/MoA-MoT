import collections

def solve_cohomology_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.

    This is based on two key results:
    1. The Künneth Formula for bounded cohomology.
    2. The structure of the bounded cohomology of Thompson's group T.
    """
    
    # Step 1: Define the dimensions of the bounded cohomology of Thompson's group T.
    # H_b^k(T; R) is an exterior algebra on a degree 3 generator.
    # We use a defaultdict to return 0 for any degree not explicitly specified.
    dim_H_T = collections.defaultdict(int)
    dim_H_T[0] = 1
    dim_H_T[3] = 1
    
    # Step 2: Set the degree n for which we want to compute the dimension.
    n = 4
    
    # Step 3: Apply the Künneth formula to compute the dimension.
    # dim H_b^n(T x T) = sum_{i=0 to n} dim H_b^i(T) * dim H_b^{n-i}(T)
    
    total_dimension = 0
    equation_terms = []
    value_terms = []
    
    for i in range(n + 1):
        j = n - i
        d_i = dim_H_T[i]
        d_j = dim_H_T[j]
        
        term_value = d_i * d_j
        total_dimension += term_value
        
        # Build strings for pretty printing the equation
        equation_terms.append(f"(dim H_b^{i}(T) * dim H_b^{j}(T))")
        value_terms.append(f"({d_i} * {d_j})")

    # Step 4: Print the results, showing the full equation.
    print(f"The dimension of the degree {n} bounded cohomology group of T x T is calculated using the Künneth formula:")
    print(f"dim H_b^{n}(T x T) = sum_{{i=0 to {n}}} dim H_b^i(T) * dim H_b^{{{n}}-i}(T)\n")

    print("The expanded formula is:")
    print(" + ".join(equation_terms))
    print("\nSubstituting the known dimensions:")
    print("= " + " + ".join(value_terms))
    print(f"= {total_dimension}")
    
    print(f"\nThe final dimension is: {total_dimension}")

if __name__ == '__main__':
    solve_cohomology_dimension()
