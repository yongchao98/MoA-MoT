import math

def solve():
    """
    Calculates the product of the first 8 terms of the sequence d_n,
    where d_n is the dimension of the B_n-fixed subspace of Kh(T(n,n); Q).
    This dimension is known to be the 2n-th Catalan number.
    """
    
    # List to store the values of d_n
    d_values = []
    
    # Calculate d_n for n from 1 to 8
    for n in range(1, 9):
        # The dimension d_n is the 2n-th Catalan number, C_{2n}.
        # The formula for the k-th Catalan number is C_k = (1/(k+1)) * C(2k, k),
        # where C(n, k) is the binomial coefficient "n choose k".
        # We calculate d_n = C_{2n} = (1/(2n+1)) * C(4n, 2n).
        # We use integer division // to avoid floating-point inaccuracies.
        k = 2 * n
        d_n = math.comb(2 * k, k) // (k + 1)
        d_values.append(d_n)
        
    # Calculate the product of all d_n values
    # Python's integers handle arbitrarily large numbers, so overflow is not an issue.
    product = 1
    for val in d_values:
        product *= val
        
    # Format the output to display the full equation
    equation_str = " * ".join(map(str, d_values))
    result_str = f"{equation_str} = {product}"
    
    # Print the final result
    print("The values of d_n for n=1 to 8 are:")
    for i, val in enumerate(d_values):
        print(f"d_{i+1} = {val}")
    
    print("\nThe product is:")
    print(result_str)

solve()