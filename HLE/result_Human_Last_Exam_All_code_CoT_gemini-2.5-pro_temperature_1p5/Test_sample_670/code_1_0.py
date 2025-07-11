def generate_tate_expression():
    """
    This function generates and prints the Tate-style program-expression graph
    for the given code snippet.
    """
    
    # Constants derived from the code snippet
    # Initial value of 'a' before the loop
    initial_a = 1
    # Initial value of 'i' in the loop
    initial_i = 1
    # Increment value for 'i' in the loop
    i_increment = 1
    
    # Print the expression part by part, outputting numbers as requested.
    # pi_1 selects the 'a' component of the result.
    print("pi_1(", end="")
    # mu is the least-fixed point operator for the loop.
    print("mu(", end="")
    # lambda defines the recurrence relation for the state (a, i).
    # 'a' and 'b' are the bound variables for the previous state.
    print("lambda(a, b). (a * b, b + ", end="")
    # The new 'i' (bound to 'b') is b + 1
    print(i_increment, end="")
    # Close the lambda body and the mu operator.
    print("))", end="")
    # Apply the mu-function to the initial state (a=1, i=1).
    print("(", end="")
    print(initial_a, end="")
    print(", ", end="")
    print(initial_i, end="")
    # Close the function application and the pi_1 projection.
    print("))")

generate_tate_expression()