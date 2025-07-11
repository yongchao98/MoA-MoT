import numpy as np
import matplotlib.pyplot as plt

def solve():
    """
    This function demonstrates that for n=2, the set S is an ellipse.
    It considers two linearly independent, non-orthogonal vectors y1, y2 in R^2,
    and generates the set of points S = {(<y1,s>^2, <y2,s>^2) | ||s||=1}.
    The resulting shape is plotted and the coefficients of its quadratic equation
    are printed.
    """
    # Define two linearly independent vectors in R^2
    y1 = np.array([1., 0.])
    y2 = np.array([2., 1.])

    # Parametrize the unit circle for s
    t = np.linspace(0, 2 * np.pi, 500)
    s_vectors = np.array([np.cos(t), np.sin(t)])

    # Compute the inner products lambda_i = <y_i, s>
    lambda1 = y1 @ s_vectors
    lambda2 = y2 @ s_vectors

    # Compute the components of the set S, x_i = lambda_i^2
    x1 = lambda1**2
    x2 = lambda2**2

    # Plot the resulting shape
    plt.figure(figsize=(6, 6))
    plt.plot(x1, x2, label=f'Shape S for y1={y1}, y2={y2}')
    plt.xlabel('$x_1 = |\langle y_1, s \\rangle|^2$')
    plt.ylabel('$x_2 = |\langle y_2, s \\rangle|^2$')
    plt.title('The shape of the set S for n=2')
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.show()

    # The general equation for the ellipse is Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0.
    # We can derive the coefficients for our specific choice of y1 and y2.
    # G = [[<y1,y1>, <y1,y2>], [<y2,y1>, <y2,y2>]]
    G = np.array([[np.dot(y1, y1), np.dot(y1, y2)],
                    [np.dot(y2, y1), np.dot(y2, y2)]])
    G_inv = np.linalg.inv(G)
    
    g11_inv, g12_inv = G_inv[0, 0], G_inv[0, 1]
    g21_inv, g22_inv = G_inv[1, 0], G_inv[1, 1]
    
    # from (g11*x1 + g22*x2 - 1)^2 = 4*(g12)^2*x1*x2
    A = g11_inv**2
    B = 2 * g11_inv * g22_inv - 4 * g12_inv**2
    C = g22_inv**2
    D = -2 * g11_inv
    E = -2 * g22_inv
    F = 1.0

    print("The shape is an ellipse described by the quadratic equation:")
    print(f"{A:.2f}*x1^2 + {B:.2f}*x1*x2 + {C:.2f}*x2^2 + {D:.2f}*x1 + {E:.2f}*x2 + {F:.2f} = 0")
    
    print("\nThe coefficients of the equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print(f"F = {F}")

solve()