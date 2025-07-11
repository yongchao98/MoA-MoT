import math

def catalan(n):
    """
    Calculates the n-th Catalan number C_n using the formula:
    C_n = (1/(n+1)) * (2n choose n)
    This is computed using integer arithmetic to avoid precision errors.
    """
    # (2n choose n) = (2n)! / (n! * n!)
    # C_n = (2n)! / ((n+1)! * n!)
    if n < 0:
        return 0
    # Using math.comb for efficient and accurate computation of combinations
    return math.comb(2 * n, n) // (n + 1)

def solve_problem():
    """
    Calculates the product of d_n for n=1 to 8.
    """
    d_values = []

    # For n=1, d_1 is the dimension of Kh(Unknot), which is 2.
    d1 = 2
    d_values.append(d1)

    # For n >= 2, d_n is the n-th Catalan number.
    for n in range(2, 9):
        dn = catalan(n)
        d_values.append(dn)

    # Calculate the final product.
    product = 1
    for d in d_values:
        product *= d

    # Format the output string to show the full equation.
    equation_str = " * ".join(map(str, d_values))
    
    print("The product is calculated as follows:")
    print(f"{equation_str} = {product}")

# Run the solver
solve_problem()