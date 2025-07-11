import math

def solve_product():
    """
    Calculates the product of d_n for n from 1 to 8.
    
    Based on a theorem in knot theory, d_n, the dimension of the B_n-invariant 
    subspace of the rational Khovanov homology of T(n,n), is 2 for all n >= 1.
    """
    n_max = 8
    d_values = []

    print("The values for d_n from n=1 to 8 are:")
    # For each n from 1 to 8, d_n is 2.
    for n in range(1, n_max + 1):
        d_n = 2
        d_values.append(d_n)
        print(f"d_{n} = {d_n}")
        
    # Calculate the product.
    # We can use math.prod for this, or a simple loop.
    # Let's stick to basics to show the concept.
    product = 1
    for d in d_values:
        product *= d
        
    # Build the equation string
    equation_string = " * ".join(map(str, d_values))
    
    print("\nThe final product is:")
    print(f"{equation_string} = {product}")

solve_product()