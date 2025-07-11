import math

def solve_for_n(n):
    """
    For a given integer n, this function calculates and prints the maximum number of
    eigenvalues greater than 2 that the matrix M(C_n) can have.
    """
    print(f"The problem is to find the maximum number of eigenvalues greater than 2 for a matrix M(C_n) of size n.")
    if not isinstance(n, int) or n < 0:
        print("Input 'n' must be a non-negative integer.")
        return

    # For n=0 and n=1, there can be no eigenvalues greater than 2.
    # The derived formula is floor((n-2)/3) + 1, which handles these cases correctly.
    if n < 2:
        result = 0
    else:
        # Using integer division for floor
        result = (n - 2) // 3 + 1
    
    # The prompt requires outputting each number in the final equation.
    print(f"\nFor n = {n}:")
    print(f"The maximum number is calculated using the formula: floor((n - 2) / 3) + 1")
    if n < 2:
         print(f"Since n < 2, the number is 0.")
    else:
        print(f"= floor(({n} - 2) / 3) + 1")
        print(f"= floor({n - 2} / 3) + 1")
        print(f"= {(n - 2) // 3} + 1")
        print(f"= {result}")
        

if __name__ == '__main__':
    # You can change this value to test for any given n.
    n_value = 14
    solve_for_n(n_value)