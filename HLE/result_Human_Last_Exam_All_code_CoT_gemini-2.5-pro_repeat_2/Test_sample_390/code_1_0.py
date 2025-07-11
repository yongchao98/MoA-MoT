import numpy as np

def get_shape_equation(y1, y2):
    """
    For n=2, calculates the coefficients of the quadratic equation
    describing the shape of the set S.
    The equation is of the form:
    A*x1^2 + C*x2^2 + B*x1*x2 + D*x1 + E*x2 + F = 0
    """
    # y1 and y2 are d-dimensional vectors
    y1 = np.array(y1)
    y2 = np.array(y2)
    
    # Check for linear independence
    # For n=2, this means they are not collinear.
    # We assume they are linearly independent as per the problem statement.
    
    # Construct the Gram matrix G
    G11 = np.dot(y1, y1)
    G12 = np.dot(y1, y2)
    G22 = np.dot(y2, y2)
    G = np.array([[G11, G12], [G12, G22]])
    
    print("Gram matrix G:")
    print(G)
    
    # Invert the Gram matrix to get Q = G^-1
    if np.linalg.det(G) == 0:
        print("Vectors are linearly dependent. Cannot proceed.")
        return

    Q = np.linalg.inv(G)
    print("\nInverse Gram matrix Q = G^-1:")
    print(Q)
    
    Q11 = Q[0, 0]
    Q12 = Q[0, 1]
    Q22 = Q[1, 1]
    
    # Coefficients of the quadratic equation for x = (x1, x2)
    # A*x1^2 + C*x2^2 + B*x1*x2 + D*x1 + E*x2 + F = 0
    A = Q11**2
    C = Q22**2
    B = 2 * Q11 * Q22 - 4 * Q12**2
    D = -2 * Q11
    E = -2 * Q22
    F = 1
    
    print("\nThe shape S is described by the quadratic equation:")
    print(f"{A:.4f}*x1^2 + {C:.4f}*x2^2 + {B:.4f}*x1*x2 + {D:.4f}*x1 + {E:.4f}*x2 + {F:.4f} = 0")
    
    # The discriminant determines the type of conic section.
    discriminant = B**2 - 4 * A * C
    
    print(f"\nThe discriminant of the quadratic form is {discriminant:.4f}.")
    if discriminant < 0:
        print("This is the equation of an ellipse (an ellipsoid in 2D).")
    elif discriminant == 0:
        print("This is the equation of a parabola (degenerates to a line/simplex in our case).")
    else:
        print("This is the equation of a hyperbola.")

# Example with a non-orthogonal basis in R^2
y1_non_ortho = [1, 0]
y2_non_ortho = [1, 1]
print("--- Analyzing non-orthogonal case ---")
get_shape_equation(y1_non_ortho, y2_non_ortho)

# Example with an orthogonal basis in R^2
y1_ortho = [2, 0]
y2_ortho = [0, 1]
print("\n\n--- Analyzing orthogonal case ---")
get_shape_equation(y1_ortho, y2_ortho)
