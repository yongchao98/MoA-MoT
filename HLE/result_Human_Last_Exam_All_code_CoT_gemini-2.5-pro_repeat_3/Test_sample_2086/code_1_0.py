import math

def solve_for_n(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n)
    based on the derived formula.

    The problem reduces to finding the maximum integer k such that n can be partitioned
    into k distinct parts, each of size at least 2. This k is given by the formula:
    k = floor((sqrt(8*n + 9) - 3) / 2)

    This function computes k for a given n and prints the intermediate steps of the calculation.

    Args:
        n: The number of objects in the category. Must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    print(f"The calculation for n = {n} is as follows:")
    
    # The final equation is k = floor((sqrt(8*n + 9) - 3) / 2)
    # We will output each number in this final equation.
    
    # Step 1: Calculate the term inside the square root: 8*n + 9
    term_in_sqrt = 8 * n + 9
    print(f"1. First, we compute 8 * n + 9: 8 * {n} + 9 = {term_in_sqrt}")
    
    # Step 2: Take the square root
    sqrt_val = math.sqrt(term_in_sqrt)
    print(f"2. Next, we take the square root: sqrt({term_in_sqrt}) = {sqrt_val}")
    
    # Step 3: Subtract 3 from the result
    numerator = sqrt_val - 3
    print(f"3. Then, we subtract 3: {sqrt_val} - 3 = {numerator}")
    
    # Step 4: Divide by 2
    division_result = numerator / 2
    print(f"4. After that, we divide by 2: {numerator} / 2 = {division_result}")
    
    # Step 5: Take the floor to get the final integer result
    max_eigenvalues = math.floor(division_result)
    print(f"5. Finally, we take the floor value: floor({division_result}) = {max_eigenvalues}")
    
    print(f"\nFor n = {n}, the maximum number of eigenvalues greater than 2 is {max_eigenvalues}.")


# You can change the value of n here to solve for a different number of objects.
n_value = 100
solve_for_n(n_value)
