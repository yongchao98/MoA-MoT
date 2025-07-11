import numpy as np

def illustrate_homotopy_on_noncompact_space():
    """
    This function illustrates a property of non-compact manifolds relevant to the problem.

    The problem's solution hinges on the manifold M being non-compact. A key property
    of non-compact manifolds is that their identity map is homotopic to a map that
    is not surjective (i.e., its image is a proper subset of the space). This
    means the space can be continuously "pushed" or "deformed" into a smaller region
    within itself.

    We use the first quadrant of the 2D plane as an example of a non-compact
    manifold M. M = {(x, y) | x > 0, y > 0}.
    
    We define a map f(p) = p + (1, 1), which translates every point.
    The image of this map is the subset U = {(x, y) | x > 1, y > 1}, which is a
    proper subset of M.
    
    This map f is homotopic to the identity map id(p) = p. The homotopy is
    H(p, t) = p + t*(1, 1), for t in [0, 1].
    H(p, 0) = p (the identity map)
    H(p, 1) = p + (1, 1) = f(p) (the non-surjective map)

    This script demonstrates this by taking sample points in M, applying the map f,
    and showing that their new positions are all within the proper subset U.
    """
    print("Illustrating a property of non-compact spaces.")
    print("Let M be the first quadrant: {(x, y) | x > 0, y > 0}.")
    print("Let f be the map f(x, y) = (x + 1, y + 1).")
    print("The image of f is the subset U = {(x, y) | x > 1, y > 1}, a proper subset of M.\n")

    # A set of sample points in our non-compact space M
    points = np.array([
        [0.5, 5.0],
        [1.0, 1.0],
        [2.3, 0.1],
        [10.0, 10.0]
    ])

    # The translation vector for our map f
    translation_vector = np.array([1.0, 1.0])

    print("Applying the map f to a set of points in M:")
    print("-" * 40)
    print(f"{'Original Point (p)':<20} | {'New Point f(p)':<20}")
    print("-" * 40)

    for p in points:
        # Apply the map f(p) = p + (1,1)
        new_p = p + translation_vector
        
        # Format for printing
        p_str = f"({p[0]:.1f}, {p[1]:.1f})"
        new_p_str = f"({new_p[0]:.1f}, {new_p[1]:.1f})"
        
        # Output each number in the final equation implicitly by showing the result
        print(f"{p_str:<20} | {new_p_str:<20}")
    
    print("-" * 40)
    print("\nAs shown, all new points (x', y') satisfy x' > 1 and y' > 1.")
    print("This demonstrates that the entire space M can be mapped into its proper subset U,")
    print("a characteristic property of non-compactness that allows the fibration to have a section.")

# Run the illustration
if __name__ == '__main__':
    illustrate_homotopy_on_noncompact_space()
