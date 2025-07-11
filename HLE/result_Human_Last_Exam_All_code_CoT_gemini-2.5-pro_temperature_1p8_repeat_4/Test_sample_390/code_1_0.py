import numpy as np
import matplotlib.pyplot as plt

def solve_and_plot_shape():
    """
    This function demonstrates that for n=2, the set S is an ellipse.
    It considers a specific case of linearly independent vectors y1, y2,
    generates points in S by sampling unit vectors s in their span,
    and plots the resulting (x1, x2) points.
    """
    # 1. Define two linearly independent vectors in R^2 (n=d=2)
    y1 = np.array([1.0, 0.0])
    y2 = np.array([1.0, 1.0])

    print("Chosen vectors:")
    print(f"y1 = {y1}")
    print(f"y2 = {y2}\n")

    # 2. The span is R^2. We sample unit vectors s from the unit circle.
    num_points = 1000
    thetas = np.linspace(0, 2 * np.pi, num_points)
    s_vectors = np.array([np.cos(thetas), np.sin(thetas)])

    # 3. For each s, calculate x = (|<y1,s>|^2, |<y2,s>|^2)
    inner_prod_1 = y1 @ s_vectors
    inner_prod_2 = y2 @ s_vectors
    x1_coords = inner_prod_1**2
    x2_coords = inner_prod_2**2

    # 4. We can derive the analytical equation for this specific case.
    # From derivation: v1 = <y1,s> = s1. x1 = s1^2
    # v2 = <y2,s> = s1+s2. x2 = (s1+s2)^2
    # s1^2+s2^2=1. Let s1=cos(theta), s2=sin(theta)
    # x1 = cos(theta)^2 = (1+cos(2*theta))/2 => 2*x1 - 1 = cos(2*theta)
    # x2 = (cos(theta)+sin(theta))^2 = 1 + 2*sin(theta)cos(theta) = 1 + sin(2*theta)
    # => x2 - 1 = sin(2*theta)
    # The equation is (2*x1-1)^2 + (x2-1)^2 = 1, which is 4*(x1-0.5)^2 + (x2-1)^2 = 1.
    # This is an ellipse centered at (0.5, 1) with semi-axes 0.5 and 1.0.

    print("The analytical equation for the shape is:")
    print("4*(x1 - 0.5)^2 + (x2 - 1.0)^2 = 1")

    # 5. Plot the results
    plt.figure(figsize=(8, 8))
    plt.plot(x1_coords, x2_coords, 'b.', markersize=4, label='Generated points from S')

    # Plot the analytical ellipse
    t = np.linspace(0, 2 * np.pi, 200)
    x1_ellipse = 0.5 + 0.5 * np.cos(t)
    x2_ellipse = 1.0 + 1.0 * np.sin(t)
    plt.plot(x1_ellipse, x2_ellipse, 'r-', linewidth=2, label='Analytical Ellipse')

    plt.title('Shape of the set S for n=2')
    plt.xlabel('$x_1 = |\langle y_1, s \\rangle|^2$')
    plt.ylabel('$x_2 = |\langle y_2, s \\rangle|^2$')
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.show()

solve_and_plot_shape()