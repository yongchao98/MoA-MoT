def solve_modulo_permutation(x, a):
    """
    Finds the permutation of a list 'a' that maximizes the final value of x
    after sequential modulo operations, and prints the calculation.

    The optimal strategy is to apply the modulo operations using the numbers
    from list 'a' in descending order.

    Args:
        x: The initial integer.
        a: A list of integers for the modulo operations.
    """
    if not a:
        print(f"The list 'a' is empty. The final x is {x}.")
        return

    # Store the original value of x to display in the equation
    original_x = x
    
    # The optimal permutation is to sort 'a' in descending order
    p = sorted(a, reverse=True)
    
    # Build the equation string step by step
    # Start with the initial x
    equation_str = f"{original_x}"
    
    # Apply the operations sequentially and update the equation string
    result = original_x
    for val in p:
        result = result % val
        equation_str = f"({equation_str} % {val})"
        
    # Print the final equation and the result
    print(f"The optimal permutation is {p}")
    print("The calculation that minimizes the difference is:")
    print(f"{equation_str} = {result}")


# --- Example Usage ---
# You can change these values to test with other inputs
initial_x = 100
number_list = [7, 12, 5, 20]

solve_modulo_permutation(initial_x, number_list)

print("\n--- Another Example ---")
initial_x_2 = 25
number_list_2 = [3, 8, 4]
solve_modulo_permutation(initial_x_2, number_list_2)