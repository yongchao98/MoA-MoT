def solve_matrix_problem(a, b):
    """
    Calculates the number of non-admissible integers k for given positive integers a and b.

    Args:
    a (int): The number of rows in the matrices. Must be a positive integer.
    b (int): The number of columns in the matrices. Must be a positive integer.
    """
    if not isinstance(a, int) or not isinstance(b, int) or a <= 0 or b <= 0:
        print("Error: a and b must be positive integers.")
        return

    print(f"Solving for a = {a}, b = {b}:")
    
    # Case 1: a = 1 or b = 1
    if a == 1 or b == 1:
        # In this case, any non-zero a x b matrix has rank 1.
        # Since all A_i are non-zero, they must all have rank 1.
        # Thus, k must be equal to ab, which is the only admissible value.
        # The non-admissible integers are 0, 1, ..., ab-1.
        ab = a * b
        result = ab
        print("Since one of the dimensions is 1, any non-zero matrix has rank 1.")
        print(f"The only admissible value for k is {ab}.")
        print(f"The number of non-admissible integers is {a} * {b} = {result}.")

    # Case 2: a > 1 and b > 1
    else:
        # For a, b > 1, an integer k is admissible if and only if it has the same parity as ab.
        ab = a * b
        
        if ab % 2 == 0:
            # If ab is even, admissible k are even. Non-admissible k are odd.
            # These are 1, 3, ..., ab-1. Their count is ab / 2.
            result = ab // 2
            print(f"The product ab = {ab} is even, so admissible k must be even.")
            print(f"The number of non-admissible integers is {ab} / 2 = {result}.")
        else: # ab is odd
            # If ab is odd, admissible k are odd. Non-admissible k are even.
            # These are 0, 2, ..., ab-1. Their count is (ab+1)/2.
            result = (ab + 1) // 2
            print(f"The product ab = {ab} is odd, so admissible k must be odd.")
            print(f"The number of non-admissible integers is ({ab} + 1) / 2 = {result}.")

# Example usage of the function.
# You can change these values to solve for different a and b.
a_val = 5
b_val = 7
solve_matrix_problem(a_val, b_val)

print("\n--- Another example ---")
a_val_2 = 4
b_val_2 = 6
solve_matrix_problem(a_val_2, b_val_2)

print("\n--- Final example ---")
a_val_3 = 1
b_val_3 = 10
solve_matrix_problem(a_val_3, b_val_3)