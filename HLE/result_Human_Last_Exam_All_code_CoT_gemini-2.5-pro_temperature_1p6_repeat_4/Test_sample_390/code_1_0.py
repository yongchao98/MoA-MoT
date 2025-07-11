import numpy as np

def solve_shape_of_S():
    """
    This function analyzes the shape of the set S for the case n=2.
    It demonstrates that for n=2, the shape is an ellipse by deriving the coefficients
    of the ellipse's equation for a specific example.

    The example uses two linearly independent, non-orthogonal vectors in R^2:
    y1 = (1, 0)
    y2 = (1, 1)

    The set S consists of points (x1, x2) where:
    x1 = |<y1, s>|^2
    x2 = |<y2, s>|^2
    and s is a unit vector in the span of y1 and y2.

    The derivation shows that x1 and x2 are related by a quadratic equation,
    which describes an ellipse. The coefficients of this equation are printed.
    """

    print("Analyzing the shape of S for a 2-dimensional example.")
    y1 = np.array([1.0, 0.0])
    y2 = np.array([1.0, 1.0])
    print(f"Using vectors y1 = {y1} and y2 = {y2}\n")

    # The shape of S is determined by an equation derived from the Gram matrix of the vectors.
    # The Gram matrix G is defined as G_ij = <y_i, y_j>.
    G = np.array([
        [np.dot(y1, y1), np.dot(y1, y2)],
        [np.dot(y2, y1), np.dot(y2, y2)]
    ])

    # The inverse of the Gram matrix, let's call it B, defines the shape.
    B = np.linalg.inv(G)
    b11, b12, b21, b22 = B.flatten()

    # For n=2, the set S is described by the points (x1, x2) which satisfy the equation:
    # (b11*x1 + b22*x2 - 1)^2 = 4 * b12^2 * x1 * x2
    # This expands to a general quadratic equation of an ellipse:
    # A*x1^2 + C*x2^2 + B_term*x1*x2 + D*x1 + E*x2 + F = 0
    
    A = b11**2
    C = b22**2
    B_term = 2 * b11 * b22 - 4 * b12**2
    D = -2 * b11
    E = -2 * b22
    F = 1.0

    print("The set S is a segment of an ellipse contained in the non-negative orthant.")
    print("The equation for this ellipse is:")
    
    equation_str = (
        f"{A:.0f}*(x1)^2 + {C:.0f}*(x2)^2 + {B_term:.0f}*(x1*x2) + "
        f"({D:.0f})*(x1) + ({E:.0f})*(x2) + {F:.0f} = 0"
    )
    
    print(equation_str)

    print("\nFor this specific example:")
    print(f"y1 = [1. 0.], y2 = [1. 1.]")
    print(f"Gram Matrix G = \n{G}")
    print(f"Inverse Gram Matrix B = \n{B}")
    print("The final equation for the ellipse is:")
    final_a = int(round(A))
    final_c = int(round(C))
    final_b_term = int(round(B_term))
    final_d = int(round(D))
    final_e = int(round(E))
    final_f = int(round(F))
    print(f"{final_a}*(x1)^2 + {final_c}*(x2)^2 + {final_b_term}*(x1*x2) + {final_d}*(x1) + {final_e}*(x2) + {final_f} = 0")


solve_shape_of_S()