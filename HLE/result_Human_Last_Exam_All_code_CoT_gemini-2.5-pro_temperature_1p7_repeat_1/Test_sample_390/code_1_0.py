import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm

def plot_S_shape():
    """
    Generates and plots points from the set S for a specific case
    of n=3 linearly independent but non-orthogonal vectors.
    """
    # Let n=3, d=3. Choose a set of linearly independent vectors y_i.
    # To make them not orthogonal, we can choose them "close" to each other.
    y1 = np.array([1, 0, 0])
    y2 = np.array([1, 1, 0])
    y3 = np.array([1, 1, 1])
    Y = np.array([y1, y2, y3]).T

    # 1. Compute the Gram matrix G = Y^T Y
    G = Y.T @ Y
    print("Gram matrix G:\n", G)

    # 2. Compute the matrix H = G^-1
    try:
        H = np.linalg.inv(G)
    except np.linalg.LinAlgError:
        print("Gram matrix is singular.")
        return
    print("\nInverse Gram matrix H:\n", H)

    # 3. Generate points on the ellipsoid u^T H u = 1.
    # We can do this by generating points x on a unit sphere,
    # and transforming them by u = H^(-1/2) x = G^(1/2) x.
    num_points = 5000
    # Generate points on a unit sphere
    x = np.random.randn(3, num_points)
    x /= np.linalg.norm(x, axis=0)

    # Compute G^(1/2) using scipy.linalg.sqrtm
    G_sqrt = sqrtm(G)
    
    # Transform points to the ellipsoid
    u = G_sqrt @ x

    # 4. Compute the points v in the set S by squaring the components of u.
    # v_i = u_i^2
    v = u**2

    # 5. Plot the resulting set S.
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(v[0, :], v[1, :], v[2, :], s=5, alpha=0.5)

    ax.set_title("Visualization of the set S for n=3 (non-orthogonal case)")
    ax.set_xlabel("$v_1 = |\\langle y_1, s \\rangle|^2$")
    ax.set_ylabel("$v_2 = |\\langle y_2, s \\rangle|^2$")
    ax.set_zlabel("$v_3 = |\\langle y_3, s \\rangle|^2$")
    ax.view_init(elev=20, azim=45)
    plt.show()

# Run the visualization
plot_S_shape()