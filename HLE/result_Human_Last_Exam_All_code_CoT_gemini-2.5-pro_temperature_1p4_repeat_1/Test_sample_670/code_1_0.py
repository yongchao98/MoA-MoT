def generate_tate_expression():
    """
    Generates and prints the Tate-style program-expression graph for 'a'
    in the given code snippet.
    """
    
    # The variable `a` is defined by a loop, which is represented by the μ (mu) operator.
    # The expression format is: μ <bound_variable> . (<initial_value>, <recurrence_expression>)

    # As per the prompt, the first bound variable is 'a'.
    bound_variable = "a"
    
    # The initial value of the variable 'a' before the loop starts.
    initial_value = 1
    
    # The recurrence expression shows how 'a' is updated in the loop.
    # It's multiplied by 'i', and the loop-carried value is represented by the bound variable 'a'.
    recurrence_expression = "a * i"
    
    # Construct the final equation string.
    final_equation = f"a = μ {bound_variable} . ({initial_value}, {recurrence_expression})"
    
    # Print the resulting equation.
    print(final_equation)

generate_tate_expression()