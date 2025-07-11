def generate_tate_expression():
    """
    This function prints the Tate-style program-expression graph for the given code.
    - π₁: Projection operator for the first element of a tuple.
    - μ(a,b): Least fixed-point operator with bound variables 'a' and 'b'.
    - φ(v1, v2): Phi function to merge values from different control paths.
    - (1, 1): The initial state of the tuple (a, i).
    - (a * b, b + 1): The update rule for the state tuple in each iteration.
    """
    # The final expression represents the value of 'a' within the loop.
    expression = "π₁(μ(a,b).(φ((1,1), (a * b, b + 1))))"
    print(expression)

generate_tate_expression()