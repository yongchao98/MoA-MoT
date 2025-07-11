def generate_tate_expression():
    """
    This function generates and prints the Tate-style program-expression
    for the specified code snippet.
    """
    # The variable 'a' is initialized to 1. This is the initial value for the loop's PHI node.
    initial_value = 1

    # Per the prompt, 'a' is the first bound variable, representing the value of the
    # variable 'a' being updated recursively in the loop.
    bound_var_a = 'a'

    # Per the prompt, 'b' is the second bound variable, representing the loop
    # counter 'i'.
    bound_var_b = 'b'

    # The recursive equation for the value of 'a' computed inside the loop is:
    # a_new = PHI(initial_value, a_new_from_previous_iteration) * i
    #
    # The μ-expression represents the least fixed-point solution to this equation.
    # We replace the recursive variable with the bound variable 'a' and the
    # loop counter 'i' with the bound variable 'b'.
    # The number `1` is the initial value in the PHI function.
    
    final_expression = f"μ {bound_var_a}.(PHI({initial_value}, {bound_var_a}) * {bound_var_b})"

    print(final_expression)

generate_tate_expression()