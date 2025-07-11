import numpy as np

def solve_max_value():
    """
    Solves for the maximum value of |b| + |c| based on the derived optimization problem.
    """
    # The problem is reduced to maximizing 2|a| + |c| subject to constraints
    # that form a parallelogram in the (a, c) plane. The maximum value must
    # occur at one of the vertices of this parallelogram.

    # Vertices of the feasible region for (a, c)
    vertices = [
        (0, 1),
        (0.5, -0.5),
        (0, -1),
        (-0.5, 0.5)
    ]

    print("The problem is to maximize |b| + |c| for the quadratic f(x) = ax^2 + bx + c,")
    print("given that |f(x)| <= 1 for all x in [-1, 1].")
    print("\nThis problem can be solved by analyzing the case where the parabola's vertex is at x = +/- 1.")
    print("This leads to an optimization problem for a and c, whose feasible region is a parallelogram.")
    print("The vertices of this parallelogram are:")
    for i, v in enumerate(vertices):
        print(f"  V{i+1}: (a, c) = {v}")

    max_val = 0
    best_params = {}

    print("\nWe evaluate the expression 2|a| + |c| (which equals |b| + |c|) at each vertex:")
    for a, c in vertices:
        # For vertex at x=1, b = -2a. For vertex at x=-1, b = 2a.
        # In both cases, |b| = 2|a|.
        b = -2 * a # We can choose b = -2a without loss of generality
        value = 2 * abs(a) + abs(c)
        print(f"  For (a, c) = ({a}, {c}), |b| + |c| = 2*|{a}| + |{c}| = {value}")
        if value > max_val:
            max_val = value
            best_params = {'a': a, 'b': b, 'c': c}

    print(f"\nThe maximum value found is {max_val}.")
    
    a = best_params['a']
    b = best_params['b']
    c = best_params['c']
    
    print("\nOne set of coefficients (a, b, c) that achieves this maximum is:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    
    print("\nThe corresponding quadratic equation is f(x) = ({a})x^2 + ({b})x + ({c}).")
    print(f"For this function, the value of |b| + |c| is |{b}| + |{c}| = {abs(b) + abs(c)}.")

solve_max_value()