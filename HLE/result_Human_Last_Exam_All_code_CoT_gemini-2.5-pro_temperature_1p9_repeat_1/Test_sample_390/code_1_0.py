import numpy as np
import matplotlib.pyplot as plt

def generate_shape_S(y_vectors, num_samples=1000):
    """
    Generates points in the set S for a given set of vectors.

    Args:
        y_vectors (list of np.array): The list of vectors [y1, y2, ...].
        num_samples (int): The number of sample points to generate.

    Returns:
        np.array: An array of points (x1, x2, ...) in the set S.
    """
    n = len(y_vectors)
    d = y_vectors[0].shape[0]

    if n != d:
        # For simplicity in this demo, we assume n=d so span({y_i}) is R^d.
        # The general logic holds for d>n as well.
        print(f"Warning: This demo works best for n=d. Currently n={n}, d={d}")
    
    # Generate random unit vectors s in the span.
    # For n=d, these are just unit vectors in R^d.
    s_vectors = np.random.randn(num_samples, d)
    s_vectors /= np.linalg.norm(s_vectors, axis=1)[:, np.newaxis]

    # Calculate the points in S
    S_points = np.zeros((num_samples, n))
    for i, s in enumerate(s_vectors):
        for j, y in enumerate(y_vectors):
            S_points[i, j] = np.abs(np.dot(y, s))**2
    
    return S_points

def main():
    """
    Main function to run the demonstration.
    """
    # Case 1: y_i are orthogonal
    # Let y1 = [2, 0] and y2 = [0, 3]. They are orthogonal but not unit vectors.
    y_ortho = [np.array([2., 0.]), np.array([0., 3.])]
    S_ortho = generate_shape_S(y_ortho)

    # For the orthogonal case, we expect x1/|y1|^2 + x2/|y2|^2 = 1
    # x1/4 + x2/9 = 1. This is a line segment (a 1-simplex).
    print("For orthogonal vectors y1=[2,0], y2=[0,3]:")
    print("The points (x1, x2) in S should satisfy the linear equation: x1/|y1|^2 + x2/|y2|^2 = 1")
    # Let's check a sample point
    sample_point = S_ortho[0]
    x1, x2 = sample_point[0], sample_point[1]
    y1_norm_sq = np.linalg.norm(y_ortho[0])**2
    y2_norm_sq = np.linalg.norm(y_ortho[1])**2
    equation_result = x1/y1_norm_sq + x2/y2_norm_sq
    print(f"Checking a random point ({x1:.4f}, {x2:.4f}):")
    print(f"{x1:.4f} / {y1_norm_sq:.1f} + {x2:.4f} / {y2_norm_sq:.1f} = {equation_result:.4f} (should be close to 1)")
    print("\n" + "="*50 + "\n")


    # Case 2: y_i are not orthogonal
    # Let y1 = [1, 0] and y2 = [1, 1]
    y_non_ortho = [np.array([1., 0.]), np.array([1., 1.])]
    S_non_ortho = generate_shape_S(y_non_ortho)
    print("For non-orthogonal vectors y1=[1,0], y2=[1,1]:")
    print("The points (x1, x2) in S should form an ellipse.")
    # The derived equation is 4*x1^2 + x2^2 - 4*x1 - 2*x2 + 1 = 0.
    x1, x2 = S_non_ortho[0][0], S_non_ortho[0][1]
    equation_result = 4*x1**2 + x2**2 - 4*x1 - 2*x2 + 1
    print(f"Checking a random point ({x1:.4f}, {x2:.4f}):")
    print(f"4*({x1:.4f})^2 + ({x2:.4f})^2 - 4*{x1:.4f} - 2*{x2:.4f} + 1 = {equation_result:.4f} (should be close to 0)")


    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot for orthogonal case
    ax1.plot(S_ortho[:, 0], S_ortho[:, 1], 'b.', markersize=2)
    ax1.set_title('Case 1: Orthogonal y_i \n The shape is a simplex')
    ax1.set_xlabel('x_1 = |<y_1, s>|^2')
    ax1.set_ylabel('x_2 = |<y_2, s>|^2')
    ax1.set_xlim(0, 4.5)
    ax1.set_ylim(0, 9.5)
    ax1.grid(True)
    ax1.set_aspect('equal', adjustable='box')

    # Plot for non-orthogonal case
    ax2.plot(S_non_ortho[:, 0], S_non_ortho[:, 1], 'r.', markersize=2)
    ax2.set_title('Case 2: Non-Orthogonal y_i \n The shape is an ellipsoid (ellipse)')
    ax2.set_xlabel('x_1 = |<y_1, s>|^2')
    ax2.set_ylabel('x_2 = |<y_2, s>|^2')
    ax2.set_xlim(-0.1, 1.1)
    ax2.set_ylim(-0.1, 2.1)
    ax2.grid(True)
    ax2.set_aspect('equal', adjustable='box')
    
    fig.suptitle('Shape of S for n=2', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == '__main__':
    main()
