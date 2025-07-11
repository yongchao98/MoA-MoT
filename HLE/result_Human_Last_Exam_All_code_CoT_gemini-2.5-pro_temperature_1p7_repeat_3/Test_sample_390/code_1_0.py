import numpy as np
import matplotlib.pyplot as plt

def solve_shape_problem():
    """
    This function visualizes the set S for a specific 2D case
    and demonstrates that the shape is an ellipse.
    """
    # Define two linearly independent vectors y1, y2
    y1 = np.array([1.0, 0.0])
    y2 = np.array([0.5, np.sqrt(3)/2.0])

    # Generate unit vectors s by varying an angle theta
    theta = np.linspace(0, 2 * np.pi, 400)
    s_vectors = np.array([np.cos(theta), np.sin(theta)]).T # A set of s vectors on the unit circle

    # Calculate the points (x1, x2) for the set S
    # x1 = |<y1, s>|^2
    # x2 = |<y2, s>|^2
    x1 = np.dot(s_vectors, y1)**2
    x2 = np.dot(s_vectors, y2)**2

    # The equation of the resulting ellipsoid can be derived.
    # Gram matrix G
    G = np.array([[np.dot(y1, y1), np.dot(y1, y2)],
                  [np.dot(y2, y1), np.dot(y2, y2)]])
    
    # Inverse of Gram matrix Q = G^{-1}
    Q = np.linalg.inv(G)
    Q11, Q12, Q22 = Q[0, 0], Q[0, 1], Q[1, 1]
    
    # Generate points for the analytically derived ellipse to compare
    # Q11^2*x^2 + Q22^2*y^2 + (2*Q11*Q22 - 4*Q12^2)*x*y - 2*Q11*x - 2*Q22*y + 1 = 0
    delta = 0.025
    x_range = np.arange(0, 1.1, delta)
    y_range = np.arange(0, 1.1, delta)
    X, Y = np.meshgrid(x_range, y_range)
    
    F = Q11**2*X**2 + Q22**2*Y**2 + (2*Q11*Q22 - 4*Q12**2)*X*Y - 2*Q11*X - 2*Q22*Y + 1
    
    # Plotting
    plt.figure(figsize=(8, 8))
    plt.plot(x1, x2, label='Set S (parametric)', color='blue', linewidth=3)
    plt.contour(X, Y, F, [0], colors='red', linestyles='dashed', linewidths=2)
    plt.title("Shape of S")
    plt.xlabel(r'$x_1 = |\langle y_1, s \rangle|^2$')
    plt.ylabel(r'$x_2 = |\langle y_2, s \rangle|^2$')
    plt.grid(True)
    plt.axis('equal')
    plt.legend(['Set S from simulation', 'Analytic Ellipse Equation'])
    plt.show()

    # We print the derived equation for the ellipsoid
    # Note: printing floating point numbers can be long.
    A = Q11**2
    B = 2*Q11*Q22 - 4*Q12**2
    C = Q22**2
    D = -2*Q11
    E = -2*Q22
    print("The derived implicit equation for the shape is of the form Ax^2 + Bxy + Cy^2 + Dx + Ey + 1 = 0")
    print(f"For y1=(1,0) and y2=(0.5, sqrt(3)/2), the equation is:")
    print(f"{A:.2f} * x1^2 + {B:.2f} * x1*x2 + {C:.2f} * x2^2 + {D:.2f} * x1 + {E:.2f} * x2 + 1.00 = 0")
    print("Since the discriminant B^2-4AC is negative, this is an ellipse.")

solve_shape_problem()