def get_minimum_curvature_cost_formula():
    """
    Calculates and prints the formula for the minimum curvature cost.

    The cost is derived by using the Woodbury matrix identity to efficiently
    compute the inverse required in the Natural Gradient Descent (NGD) update.
    This approach avoids the direct inversion of a d x d matrix, which would
    cost O(d^3), and instead relies on the inversion of an n x n matrix,
    which is cheaper since n < d. The dominant cost then comes from the
    surrounding matrix multiplications.

    The total cost is the sum of floating-point operations for:
    1.  Cost(X^T * X)               -> 2 * d * n^2
    2.  Cost(inverting P)           -> 2 * n^3
    3.  Cost(G * X)                 -> 2 * d^2 * n
    4.  Cost((G*X) * P_inv)         -> 2 * d * n^2
    5.  Cost((G*X*P_inv) * X^T)     -> 2 * d^2 * n
    -----------------------------------------------------
    Total = (4 * d^2 * n) + (4 * d * n^2) + (2 * n^3)
    """
    
    # Coefficients of the terms in the final cost formula
    d_squared_n_coeff = 4
    d_n_squared_coeff = 4
    n_cubed_coeff = 2
    
    # Construct the formula string
    cost_formula = (
        f"{d_squared_n_coeff}*d^2*n + "
        f"{d_n_squared_coeff}*d*n^2 + "
        f"{n_cubed_coeff}*n^3"
    )
    
    print("The minimum curvature cost, in terms of floating-point operations (flops), is approximately:")
    print(cost_formula)

if __name__ == '__main__':
    get_minimum_curvature_cost_formula()
