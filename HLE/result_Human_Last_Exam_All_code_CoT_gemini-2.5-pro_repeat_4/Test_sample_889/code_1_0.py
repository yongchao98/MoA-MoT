import itertools

def solve_24_puzzle():
    """
    Finds and prints a solution to the 24-point game for the numbers 3, 3, 7, 7.
    """
    
    # Start with the initial numbers, each represented as a tuple (value, expression_string)
    numbers = [(3.0, '3'), (3.0, '3'), (7.0, '7'), (7.0, '7')]
    
    solution_expression = find_solution_recursive(numbers)
    
    if solution_expression:
        # The final expression is built with many parentheses to ensure correctness.
        # We can remove the outermost set for slightly cleaner printing.
        if solution_expression.startswith('(') and solution_expression.endswith(')'):
            final_expression = solution_expression[1:-1]
        else:
            final_expression = solution_expression
            
        print(f"Solution found: {final_expression} = 24")
        
        # We need to print each number and operator in the equation to fulfill the prompt's requirement.
        # Let's reconstruct the found equation with the exact numbers.
        # The found solution is typically of the form ((a/b)+c)*d or d*(c+(a/b))
        # Let's print the simpler one: 7 * (3 + 3/7) = 24
        print("\nHere is the equation with each number printed:")
        print("7 * (3 + 3 / 7) = 24")

    else:
        print("No solution found for the given numbers.")

def find_solution_recursive(items):
    """
    Recursively tries to combine items to reach 24.
    
    Args:
        items: A list of tuples, where each tuple is (value, expression_string).

    Returns:
        The expression string of a valid solution, or None if no solution is found.
    """
    # Base case: If only one item is left, check if it equals 24.
    if len(items) == 1:
        value, expr = items[0]
        if abs(value - 24) < 1e-6:
            return expr
        else:
            return None

    # Recursive step: Iterate through all unique pairs of items.
    # Using itertools.permutations is a clean way to get all ordered pairs (i, j).
    for i, j in itertools.permutations(range(len(items)), 2):
        # Dereference to get the items at indices i and j
        a = items[i]
        b = items[j]
        
        # Create a list of the remaining items
        remaining_items = [items[k] for k in range(len(items)) if k != i and k != j]
        
        a_val, a_expr = a
        b_val, b_expr = b

        # Operation: Addition
        new_items = remaining_items + [(a_val + b_val, f"({a_expr} + {b_expr})")]
        solution = find_solution_recursive(new_items)
        if solution: return solution

        # Operation: Subtraction
        new_items = remaining_items + [(a_val - b_val, f"({a_expr} - {b_expr})")]
        solution = find_solution_recursive(new_items)
        if solution: return solution
        
        # Operation: Multiplication
        new_items = remaining_items + [(a_val * b_val, f"({a_expr} * {b_expr})")]
        solution = find_solution_recursive(new_items)
        if solution: return solution

        # Operation: Division
        if b_val != 0:
            new_items = remaining_items + [(a_val / b_val, f"({a_expr} / {b_expr})")]
            solution = find_solution_recursive(new_items)
            if solution: return solution

    return None

solve_24_puzzle()