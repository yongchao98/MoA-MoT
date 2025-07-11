def generate_tate_expression_graph():
    """
    This function generates and prints the Tate-style program-expression graph for the given code snippet.

    The code is:
    int a = 1;
    for (int i = 1; i < n; i++) { a *= i; }

    1. The value of the loop counter 'i' starts at 1 and increments by 1 each time.
       Its program-expression graph is `μ b . (1 | b + 1)`.
       We use 'b' as the second bound variable as requested.

    2. The value of 'a' starts at 1 and in each iteration is updated to `a * i`.
       Its program-expression graph depends on its own previous value (the bound variable 'a')
       and the value of 'i'. The graph is `μ a . (1 | a * [expression_for_i])`.
       We use 'a' as the first bound variable.

    3. By substituting the expression for 'i' into the expression for 'a', we get the
       complete expression graph for the value of 'a' in the loop.
    """

    # Expression for the loop counter 'i'
    # Initial value is 1, recurrence is b + 1. Bound variable is b.
    i_expression = "μ b . (1 | b + 1)"

    # Expression for the variable 'a'
    # Initial value is 1, recurrence is a * (value of i). Bound variable is a.
    # We embed the expression for 'i' into this one.
    a_expression = f"μ a . (1 | a * ({i_expression}))"

    print(a_expression)

if __name__ == "__main__":
    generate_tate_expression_graph()