def solve():
    """
    This function generates the Tate-style program-expression graph for the given code snippet.

    The C code is:
    int a = 1;
    for (int i = 1; i < n; i++) { a *= i; }

    Analysis:
    1. Initial value: The variable 'a' is initialized to 1.
    2. Update rule: Inside the loop, 'a' is updated by 'a * i'.
    3. Bound variables: The problem specifies using 'a' for the first bound variable 
       (the loop-carried value of 'a') and 'b' for the second (the value of 'i').

    Tate-style expression format for a loop is: mu <var> . (<initial>, <update_expr>)
    - <var> is 'a'.
    - <initial> is 1.
    - <update_expr> is 'a * b'.
    """

    # Define the components of the Tate expression based on the analysis
    initial_value = 1
    first_bound_var = "a"
    # The update expression `a * i` becomes `a * b` using the specified bound variables
    update_expression = "a * b"

    # Construct the final expression string, ensuring the number is outputted directly.
    # The final expression represents the value of the variable 'a' defined by the loop.
    final_expression = f"μ {first_bound_var} . ({initial_value}, {update_expression})"

    print(final_expression)

solve()