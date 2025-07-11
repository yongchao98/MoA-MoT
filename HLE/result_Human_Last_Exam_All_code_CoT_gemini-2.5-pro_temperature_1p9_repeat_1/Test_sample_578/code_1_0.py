import math

def solve():
    """
    This function calculates the product of the first 8 terms of the sequence d_n,
    where d_n is the dimension of the subspace of Kh(T(n,n); Q) fixed by B_n.
    Based on established results, d_n = n + 1.
    The problem asks for the product of d_n for n from 1 to 8.
    """
    
    # Calculate d_n for n=1 to 8
    d_values = [n + 1 for n in range(1, 9)]
    
    # Calculate the product of d_n values
    product = 1
    for d in d_values:
        product *= d
        
    # The user wants each number in the final equation printed.
    # Create the string for the equation.
    equation_str = " * ".join(map(str, d_values))
    
    # Print the equation and the final result.
    print(f"The values of d_n for n=1 to 8 are: {d_values}")
    print(f"The product is calculated as:")
    print(f"{equation_str} = {product}")

solve()