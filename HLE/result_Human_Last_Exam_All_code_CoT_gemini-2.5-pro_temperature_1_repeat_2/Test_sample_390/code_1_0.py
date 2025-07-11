import numpy as np

def demonstrate_shapes():
    """
    Demonstrates that the shape of the set S depends on the choice of vectors y_i.
    """
    # Case 1: Orthogonal vectors y_i
    print("Case 1: Orthogonal vectors y1 = [1, 0], y2 = [0, 2]")
    y1_ortho = np.array([1, 0])
    y2_ortho = np.array([0, 2])
    
    print("The points (x1, x2) in S satisfy the linear equation: a*x1 + b*x2 = c")
    # Equation: x1 + x2/4 = 1
    a1, b1, c1 = 1, 0.25, 1
    print(f"The equation is: {a1}*x1 + {b1}*x2 = {c1}")
    print("This is the equation of a line segment (a 1-simplex).\n")
    
    print("Sample points from S for the orthogonal case:")
    points_ortho = []
    for theta in np.linspace(0, np.pi/2, 5): # a quarter is enough to see the segment
        s = np.array([np.cos(theta), np.sin(theta)])
        x1 = np.dot(y1_ortho, s)**2
        x2 = np.dot(y2_ortho, s)**2
        points_ortho.append((x1, x2))
        print(f"  theta={theta:.2f}, (x1, x2) = ({x1:.4f}, {x2:.4f})")
        
    print("-" * 30)

    # Case 2: Non-orthogonal vectors y_i
    print("Case 2: Non-orthogonal vectors y1 = [1, 0], y2 = [1, 1]")
    y1_non_ortho = np.array([1, 0])
    y2_non_ortho = np.array([1, 1])

    print("The points (x1, x2) in S satisfy the quadratic equation: a*(x1-h)^2 + b*(x2-k)^2 = r^2")
    # Equation: 4*(x1 - 0.5)^2 + (x2 - 1)^2 = 1
    a2, h2, b2, k2, r2_sq = 4, 0.5, 1, 1, 1
    print(f"The equation is: {a2}*(x1 - {h2})^2 + {b2}*(x2 - {k2})^2 = {r2_sq}")
    print("This is the equation of an ellipse.\n")

    print("Sample points from S for the non-orthogonal case:")
    points_non_ortho = []
    for theta in np.linspace(0, 2 * np.pi, 9)[:-1]: # a full circle
        s = np.array([np.cos(theta), np.sin(theta)])
        x1 = np.dot(y1_non_ortho, s)**2
        x2 = np.dot(y2_non_ortho, s)**2
        points_non_ortho.append((x1, x2))
        print(f"  theta={theta:.2f}, (x1, x2) = ({x1:.4f}, {x2:.4f})")

demonstrate_shapes()