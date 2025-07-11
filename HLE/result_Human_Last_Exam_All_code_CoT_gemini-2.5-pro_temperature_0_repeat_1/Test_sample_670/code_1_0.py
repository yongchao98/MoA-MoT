def generate_tate_expression():
    """
    Generates and prints the Tate-style program-expression for the given loop.
    """
    # The state is (a, i), represented by bound variables (a, b).
    # Initial state is a=1, i=1, so (1, 1).
    initial_state = "(1, 1)"
    
    # The loop condition is i < n, which is b < n.
    condition = "b < n"
    
    # The state update is a = a * i and i = i + 1, so (a * b, b + 1).
    update_expression = "(a * b, b + 1)"
    
    # The recursive expression for the loop body.
    loop_body = f"if {condition} then {update_expression} else (a, b)"
    
    # The mu-expression for the final state.
    mu_expression = f"(μ a,b. ({loop_body}))"
    
    # The final state is the mu-expression applied to the initial state.
    final_state = f"{mu_expression}{initial_state}"
    
    # We want the expression for 'a', which is the first component of the state.
    # We use π₁ to denote projection onto the first component.
    final_expression_for_a = f"π₁({final_state})"
    
    print(final_expression_for_a)

generate_tate_expression()