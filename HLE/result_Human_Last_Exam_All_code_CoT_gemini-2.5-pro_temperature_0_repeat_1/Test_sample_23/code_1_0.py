import math

def count_non_admissible_integers(a, b):
    """
    Calculates and explains the number of non-admissible integers k for given a and b.

    Args:
        a (int): The number of rows in the matrices, a positive integer.
        b (int): The number of columns in the matrices, a positive integer.
    """
    print(f"Solving for a = {a} and b = {b}:")
    print("=" * 30)

    # Case 1: One of the dimensions is 1.
    if a == 1 or b == 1:
        print("Analysis: Since one dimension is 1, any non-zero matrix is a vector and has rank 1.")
        print("All basis matrices A_i must be non-zero, so all must have rank 1.")
        ab = a * b
        print(f"This means k must be equal to a*b = {ab}.")
        print(f"The only admissible integer is k = {ab}.")
        print(f"The non-admissible integers are all other values in the range [0, {ab}], which are {{0, 1, ..., {ab-1}}}.")
        
        num_not_admissible = ab
        
        print("\nFinal Calculation:")
        print(f"Number of non-admissible integers = {a} * {b} = {num_not_admissible}")

    # Case 2: Both dimensions are greater than 1.
    else:
        print("Analysis: Since a > 1 and b > 1, we can use the parity theorem.")
        print("The number of rank-1 matrices, k, in any orthogonal basis must have the same parity.")
        ab = a * b
        print(f"The standard basis (E_ij) has {ab} rank-1 matrices. So, k must have the same parity as {ab}.")
        
        if ab % 2 == 0:
            print(f"a*b = {ab} is even. Admissible k must be even.")
            print(f"The non-admissible integers are the odd numbers in the range [0, {ab}].")
            print("\nFinal Calculation:")
            result = ab // 2
            print(f"Number of non-admissible integers = {ab} / 2 = {result}")
        else:
            print(f"a*b = {ab} is odd. Admissible k must be odd.")
            print(f"The non-admissible integers are the even numbers in the range [0, {ab}].")
            print("\nFinal Calculation:")
            result = (ab + 1) // 2
            print(f"Number of non-admissible integers = ({ab} + 1) / 2 = {result}")

# --- Main execution ---
# Let's solve for a specific case, e.g., a=10, b=13
a_val = 10
b_val = 13
count_non_admissible_integers(a_val, b_val)