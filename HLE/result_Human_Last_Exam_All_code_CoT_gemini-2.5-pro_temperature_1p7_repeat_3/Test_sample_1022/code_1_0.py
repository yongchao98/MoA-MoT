import math

def compute_dimension():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T
    using the Künneth formula and known results for Thompson's group T.
    """
    
    # Step 1: Define the known dimensions for H_b^n(T, R).
    # We use math.inf for infinity. Dimensions for n>=4 are 0.
    dims_H_b_T = {
        0: 1,
        1: 0,
        2: math.inf,
        3: 1,
    }

    # Step 2: Set the degree for the calculation.
    degree = 4
    
    term_values = []
    calculation_steps = []

    # Step 3: Loop through all pairs (p, q) such that p + q = degree.
    for p in range(degree + 1):
        q = degree - p
        
        # Get dimensions for p and q, defaulting to 0 for degrees >= 4.
        dim_p = dims_H_b_T.get(p, 0)
        dim_q = dims_H_b_T.get(q, 0)
        
        # In the context of dimensions, 0 * infinity = 0.
        if dim_p == 0 or dim_q == 0:
            term_product = 0
        else:
            term_product = dim_p * dim_q
        
        term_values.append(term_product)
        
        # Prepare strings for clear output.
        dim_p_str = "∞" if dim_p == math.inf else str(dim_p)
        dim_q_str = "∞" if dim_q == math.inf else str(dim_q)
        product_str = "∞" if term_product == math.inf else str(int(term_product))
        calculation_steps.append(f"p={p}, q={q}: dim(H_b^{p})*dim(H_b^{q}) = {dim_p_str} * {dim_q_str} = {product_str}")

    # Step 4: Sum the products to get the total dimension.
    total_dimension = sum(term_values)
    final_result_str = "∞" if total_dimension == math.inf else str(int(total_dimension))
    
    # Step 5: Print the detailed calculation.
    print("Computing dim H_b^4(T x T, R) using the Künneth formula:")
    print("dim = sum_{p+q=4} [dim H_b^p(T, R) * dim H_b^q(T, R)]\n")
    
    print("Calculation for each term (p,q):")
    for step in calculation_steps:
        print(step)
    
    print("\nFinal Equation:")
    # This line prints each number in the final equation as requested.
    sum_equation = " + ".join(["∞" if v == math.inf else str(int(v)) for v in term_values])
    print(f"dim H_b^4(T x T, R) = {sum_equation}")
    print(f"                     = {final_result_str}")

if __name__ == "__main__":
    compute_dimension()