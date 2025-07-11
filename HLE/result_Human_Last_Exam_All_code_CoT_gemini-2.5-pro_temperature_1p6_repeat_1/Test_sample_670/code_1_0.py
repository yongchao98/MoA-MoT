import sys

def solve():
    """
    This function generates and prints the Tate-style program-expression graph
    for the given code snippet.

    The C++ code is:
    int a = 1;
    for (int i = 1; i < n; i++) { a *= i; }

    1. The loop counter 'i' starts at 1 and has the recurrence i -> i + 1.
       Using 'b' as the bound variable, its expression is mu b.(1 + b).

    2. The variable 'a' starts at 1 and has the recurrence a -> a * i.
       Its recurrence depends on the value of 'i' at each step of the loop.
       Using 'a' as the bound variable, we express this by nesting the
       expression for 'i' inside a's mu-expression.

    3. The resulting expression for 'a' is mu a.(a * (mu b.(1 + b))).
    """

    # The expression for the loop counter 'i', using 'b' as the bound variable.
    # The recurrence is i -> i + 1. The number is 1.
    inner_expression = "μ b.(1 + b)"

    # The expression for 'a', using 'a' as the bound variable.
    # The recurrence is a -> a * i.
    # We embed the expression for 'i' to represent the dependency.
    final_expression = f"μ a.(a * ({inner_expression}))"

    print(final_expression)

solve()