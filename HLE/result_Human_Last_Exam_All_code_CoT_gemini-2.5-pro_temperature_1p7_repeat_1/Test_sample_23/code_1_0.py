def solve_and_explain(a, b):
    """
    Calculates the number of non-admissible integers k for given positive integers a and b.

    Args:
        a: The number of rows in the matrices.
        b: The number of columns in the matrices.
    """
    print(f"For a = {a} and b = {b}:")

    # Case 1: One of the dimensions is 1.
    if min(a, b) == 1:
        result = a * b
        print("Since one of the dimensions is 1, any non-zero matrix has rank 1.")
        print("This forces any basis of non-zero matrices to consist entirely of rank-1 matrices.")
        print("Thus, only k=a*b is admissible.")
        print("The number of non-admissible integers in {0, 1, ..., a*b} is a*b.")
        print(f"The calculation is: {a} * {b} = {result}")

    # Case 2: Both dimensions are 2 or greater.
    else:  # a, b >= 2
        result = a + b - 2
        print("Since both dimensions are 2 or greater, we are in the more complex case.")
        print("Based on known results and the geometry of rank-1 matrices, the number of non-admissible integers is a + b - 2.")
        print(f"The calculation is: {a} + {b} - 2 = {result}")

# --- Main execution ---
# We demonstrate the function with a few examples.

# Example where a,b >= 2
solve_and_explain(3, 4)
print("-" * 20)

# Example where min(a,b) = 1
solve_and_explain(1, 5)
print("-" * 20)

# Example for the well-understood a=2, b=2 case
solve_and_explain(2, 2)