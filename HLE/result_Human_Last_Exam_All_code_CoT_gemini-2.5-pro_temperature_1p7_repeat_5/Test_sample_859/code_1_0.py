import math

def solve_graph_problem():
    """
    This function explains and calculates the minimal number of new edges required.
    Since the degree 'd' is symbolic, the function outputs the answer as a formula in terms of 'd'.
    """

    # Representing 'd' symbolically for the output expressions
    d_str = 'd'
    
    # The numbers used in the formulas
    one = 1
    two = 2
    three = 3
    four = 4
    six = 6
    
    print(f"The minimal number of new edges, N, depends on the even integer {d_str}, where {d_str} >= 2.")

    # Case for d = 2
    print(f"\nIf {d_str} == {two}:")
    num_leaves_d2 = 6
    result_d2 = 3
    print(f"The maximum number of leaves that can be formed in G' is {num_leaves_d2}.")
    print(f"The number of new edges to add is N = ceil({num_leaves_d2} / {two}) = {result_d2}.")

    # Case for d >= 4
    print(f"\nIf {d_str} >= {four} (and is an even integer):")
    print(f"The maximum number of leaves that can be formed in G' is {three}*{d_str} + {two}.")
    print(f"The number of new edges to add is N = ceil(({three}*{d_str} + {two}) / {two}) = ({three}*{d_str} / {two}) + {one}.")

solve_graph_problem()
# The final answer is a piecewise formula depending on d.
final_answer_formula = "3 if d == 2 else (3*d/2 + 1)"