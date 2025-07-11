import numpy as np

def get_shape_equation(y1, y2):
    """
    For n=2, derives the coefficients of the quadratic equation that defines the shape S.
    The equation is Ax_1^2 + Bx_1x_2 + Cx_2^2 + Dx_1 + Ex_2 + F = 0.
    """
    # y1 and y2 are given as numpy arrays
    y1 = np.array(y1, dtype=float)
    y2 = np.array(y2, dtype=float)
    
    # Gram matrix G
    G = np.array([
        [np.dot(y1, y1), np.dot(y1, y2)],
        [np.dot(y2, y1), np.dot(y2, y2)]
    ])
    
    # Check for linear independence
    if np.linalg.det(G) < 1e-9:
        print("Vectors are not linearly independent.")
        return

    # H is the inverse of the Gram matrix
    H = np.linalg.inv(G)
    h11, h12, h22 = H[0, 0], H[0, 1], H[1, 1]

    # Coefficients of the ellipse equation:
    # A*x1^2 + C*x2^2 + B*x1*x2 + D*x1 + E*x2 + F = 0
    A = h11**2
    C = h22**2
    B = 2 * h11 * h22 - 4 * h12**2
    D = -2 * h11
    E = -2 * h22
    F = 1.0

    print(f"For y1={y1.tolist()}, y2={y2.tolist()}:")
    print("The shape S is an ellipse defined by the quadratic equation:")
    print(f"{A:.4f}*x1^2 + {B:.4f}*x1*x2 + {C:.4f}*x2^2 + {D:.4f}*x1 + {E:.4f}*x2 + {F:.4f} = 0")

# Case 1: Orthogonal vectors (leads to a simplex, a degenerate ellipse)
y1_ortho = [2.0, 0.0]
y2_ortho = [0.0, 3.0]
get_shape_equation(y1_ortho, y2_ortho)
# For orthogonal vectors, G and H are diagonal, so h12=0. The equation simplifies to
# h11^2*x1^2 + h22^2*x2^2 + 2*h11*h22*x1*x2 - 2*h11*x1 - 2*h22*x2 + 1 = 0
# which is (h11*x1 + h22*x2 - 1)^2 = 0, so h11*x1 + h22*x2 = 1.
# Since H=G^-1, h11=1/||y1||^2 and h22=1/||y2||^2, so x1/||y1||^2 + x2/||y2||^2 = 1.
# This is a 1-simplex (a line segment).
print("-" * 20)

# Case 2: Non-orthogonal vectors (generic case)
# using y2 rotated by 60 degrees from y1
angle = np.pi / 3
y1_non_ortho = [1.0, 0.0]
y2_non_ortho = [np.cos(angle), np.sin(angle)]
get_shape_equation(y1_non_ortho, y2_non_ortho)
