import numpy as np
import matplotlib.pyplot as plt

def analyze_shape(y_vectors, title):
    """
    Generates and plots the set S for a given set of 2 vectors in 2D.
    It also computes and prints the equation that describes the shape.
    """
    y1 = np.array(y_vectors[0])
    y2 = np.array(y_vectors[1])

    # The span of two linearly independent vectors in R^2 is R^2.
    # We generate random unit vectors s in R^2.
    thetas = np.linspace(0, 2 * np.pi, 1000)
    s_vectors = np.array([np.cos(thetas), np.sin(thetas)])

    # Calculate the points in S
    x1_coords = np.square(np.dot(y1, s_vectors))
    x2_coords = np.square(np.dot(y2, s_vectors))

    # Plot the shape of S
    plt.figure(figsize=(8, 8))
    plt.plot(x1_coords, x2_coords, label='Set S')
    plt.title(title)
    plt.xlabel(r'$x_1 = |\langle y_1, s \rangle|^2$')
    plt.ylabel(r'$x_2 = |\langle y_2, s \rangle|^2$')
    plt.grid(True)
    plt.axis('equal')
    
    # Calculate Gram matrix and its inverse
    G = np.array([
        [np.dot(y1, y1), np.dot(y1, y2)],
        [np.dot(y2, y1), np.dot(y2, y2)]
    ])
    
    print(f"--- {title} ---")
    if np.abs(G[0, 1]) < 1e-9: # Orthogonal case
        # Equation is x1/||y1||^2 + x2/||y2||^2 = 1
        c1 = 1 / G[0, 0]
        c2 = 1 / G[1, 1]
        print("Shape is a simplex (line segment).")
        print(f"Equation: ({c1:.4f})*x1 + ({c2:.4f})*x2 = 1")
        x_theory = np.linspace(0, G[0,0], 100)
        y_theory = G[1,1] * (1 - x_theory/G[0,0])
        plt.plot(x_theory, y_theory, 'r--', label='Theoretical Simplex')

    else: # Non-orthogonal case
        # Equation from (g11*x1 + g22*x2 - 1)^2 = 4*g12^2*x1*x2
        # A*x1^2 + B*x2^2 + C*x1*x2 + D*x1 + E*x2 + F = 0
        G_inv = np.linalg.inv(G)
        g11, g12, g22 = G_inv[0, 0], G_inv[0, 1], G_inv[1, 1]
        A = g11**2
        B = g22**2
        C = 2*g11*g22 - 4*g12**2
        D = -2*g11
        E = -2*g22
        F = 1
        print("Shape is an ellipse.")
        print("Equation of the ellipse bounding the shape:")
        print(f"({A:.4f})x1^2 + ({B:.4f})x2^2 + ({C:.4f})x1*x2 + ({D:.4f})x1 + ({E:.4f})x2 + {F:.4f} = 0")
        
        # For plotting the theoretical curve
        delta = 0.005
        x_range = np.arange(np.min(x1_coords)-0.1, np.max(x1_coords)+0.1, delta)
        y_range = np.arange(np.min(x2_coords)-0.1, np.max(x2_coords)+0.1, delta)
        X, Y = np.meshgrid(x_range,y_range)
        Z = A*X**2 + B*Y**2 + C*X*Y + D*X + E*Y + F
        plt.contour(X, Y, Z, [0], colors='r', linestyles='dashed')
        plt.legend()
        
    plt.show()


# Case 1: Orthogonal vectors
y_orth = [(2.0, 0.0), (0.0, 3.0)]
analyze_shape(y_orth, 'Case 1: Orthogonal Vectors')

# Case 2: Non-orthogonal vectors
y_non_orth = [(1.0, 0.0), (1.0, 1.0)]
analyze_shape(y_non_orth, 'Case 2: Non-Orthogonal Vectors')
