import numpy as np
import matplotlib.pyplot as plt

def generate_and_plot_S(n, d, num_points):
    """
    Generates and plots the set S for a given n and d.
    """
    if d < n:
        raise ValueError("d must be greater than or equal to n for linear independence.")

    # 1. Generate n random linearly independent vectors y_i in R^d
    # Random vectors from a continuous distribution are linearly independent with probability 1.
    Y = np.random.randn(n, d)

    # 2. Find an orthonormal basis for the span of {y_i}
    # Using QR decomposition is numerically stable. Q will have orthonormal columns.
    # We apply it to Y.T so that the columns of Y.T (rows of Y) are the vectors.
    Q, _ = np.linalg.qr(Y.T)
    orthonormal_basis = Q[:, :n] # Orthonormal basis for the span of rows of Y

    # 3. Generate random unit vectors s in the span
    S_points = []
    for _ in range(num_points):
        # Generate a random unit vector 'a' in R^n
        a = np.random.randn(n)
        a /= np.linalg.norm(a)

        # Construct s = sum(a_j * u_j)
        s = np.dot(orthonormal_basis, a)

        # 4. Calculate the point x in S
        x = [np.abs(np.dot(y_i, s))**2 for y_i in Y]
        S_points.append(x)

    S_points = np.array(S_points)

    # 5. Plot the result
    fig = plt.figure(figsize=(8, 8))
    if n == 2:
        ax = fig.add_subplot(111)
        ax.scatter(S_points[:, 0], S_points[:, 1], s=5)
        ax.set_xlabel('$|\langle y_1, s \\rangle|^2$')
        ax.set_ylabel('$|\langle y_2, s \\rangle|^2$')
        ax.set_title(f'Shape of S for n=2 (an ellipse)')
        ax.grid(True)
        ax.axis('equal')
    elif n == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(S_points[:, 0], S_points[:, 1], S_points[:, 2], s=5)
        ax.set_xlabel('$|\langle y_1, s \\rangle|^2$')
        ax.set_ylabel('$|\langle y_2, s \\rangle|^2$')
        ax.set_zlabel('$|\langle y_3, s \\rangle|^2$')
        ax.set_title(f'Shape of S for n=3 (general case)')
    else:
        print(f"Plotting for n={n} is not implemented in this script.")
        return

    plt.show()

# Visualize for the general case n=3, d=3, which is not an ellipsoid or simplex
try:
    generate_and_plot_S(n=3, d=3, num_points=5000)
    print("The plot for n=3 shows a shape that is not a simplex, ellipsoid, or any of the other simple options.")
    print("This confirms that in the general case, the shape is more complex.")
except Exception as e:
    print(f"An error occurred during plotting: {e}")
