import numpy as np
from scipy.optimize import linprog
import math

def solve_max_ghtw():
    """
    Calculates the maximum generalized hypertreewidth for a hypergraph with 3 edges.
    
    The problem is solved by finding the maximum fractional hypertreewidth (fhw)
    and then taking the ceiling of that value. The fhw is found by solving the
    dual linear programming problem, which is equivalent to finding the
    fractional matching number alpha*.

    We demonstrate this with the worst-case hypergraph structure:
    E = {{a, d}, {b, d}, {c, d}}
    Vertices are represented as an array index [a, b, c, d].
    
    The LP is:
    Maximize: y_a + y_b + y_c + y_d
    Subject to:
        y_a + y_d <= 1
        y_b + y_d <= 1
        y_c + y_d <= 1
        y_i >= 0
    
    `linprog` minimizes, so we minimize the negative of the objective function.
    """
    
    # Objective function to minimize: -y_a - y_b - y_c - y_d
    c = [-1, -1, -1, -1]

    # The inequality constraints matrix A_ub (LHS of Ax <= b)
    # y_a + y_d <= 1  -> [1, 0, 0, 1]
    # y_b + y_d <= 1  -> [0, 1, 0, 1]
    # y_c + y_d <= 1  -> [0, 0, 1, 1]
    A_ub = np.array([
        [1, 0, 0, 1],
        [0, 1, 0, 1],
        [0, 0, 1, 1]
    ])

    # The inequality constraints vector b_ub (RHS of Ax <= b)
    b_ub = np.array([1, 1, 1])

    # Bounds for each variable (y_i >= 0)
    bounds = (0, None)

    # Solve the linear program
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

    if result.success:
        # The result.fun is the minimum value of the objective function (-z).
        # The maximum value (alpha*) is -result.fun.
        max_alpha_star = -result.fun
        
        # The max fhw is equal to the max alpha*
        max_fhw = max_alpha_star
        
        # The max ghtw is the ceiling of the max fhw
        max_ghtw = math.ceil(max_fhw)
        
        print(f"The maximum fractional matching number (α*) is: {max_alpha_star}")
        print(f"This means the maximum fractional hypertreewidth (fhw) is: {max_fhw}")
        print(f"The maximum generalized hypertreewidth (ghtw) is ceil(fhw), which is: ceil({max_fhw}) = {max_ghtw}")

    else:
        print("The linear programming optimization failed.")
        print(result.message)

solve_max_ghtw()