import numpy as np

def solve():
    """
    Finds the largest possible length of a line segment on the unit sphere
    of a space with a metrically convex unit ball.
    """
    print("We will demonstrate the solution in the space R^2 with the L_infinity norm.")
    print("This space's unit ball is metrically convex.")
    print("The unit sphere is the boundary of the square [-1, 1] x [-1, 1].\n")

    # 1. Define two points on the unit sphere.
    # We choose two vertices of the square.
    u = np.array([1, 1])
    v = np.array([1, -1])

    # Check that they are on the unit sphere using the L_infinity norm.
    norm_u = np.linalg.norm(u, np.inf)
    norm_v = np.linalg.norm(v, np.inf)

    print(f"Chosen point u = ({u[0]}, {u[1]}) with L_infinity norm: {norm_u}")
    print(f"Chosen point v = ({v[0]}, {v[1]}) with L_infinity norm: {norm_v}\n")

    # 2. Verify that the line segment between u and v lies on the unit sphere.
    # A point on the segment is x_t = (1-t)*u + t*v for t in [0, 1].
    # Analytically, x_t = (1, 1-2t).
    # The norm is ||x_t||_inf = max(|1|, |1-2t|).
    # For t in [0, 1], |1-2t| <= 1, so the norm is always 1.
    print("Verifying that the segment lies on the unit sphere:")
    is_on_sphere = True
    for t in np.linspace(0, 1, 5):
        x_t = (1 - t) * u + t * v
        norm_x_t = np.linalg.norm(x_t, np.inf)
        print(f"  For t = {t:.2f}, point x_t = ({x_t[0]:.2f}, {x_t[1]:.2f}), norm = {norm_x_t:.2f}")
        if not np.isclose(norm_x_t, 1.0):
            is_on_sphere = False
    
    if is_on_sphere:
        print("Conclusion: The entire segment lies on the unit sphere.\n")
    else:
        print("Error: The segment does not lie on the unit sphere.\n")

    # 3. Calculate the length of the segment.
    # The length is the L_infinity norm of the difference vector.
    diff_vec = u - v
    length = np.linalg.norm(diff_vec, np.inf)

    print("The length of the line segment is calculated as ||u - v||_inf.")
    # As requested, outputting the numbers in the final equation.
    print(f"Length = ||({u[0]}, {u[1]}) - ({v[0]}, {v[1]})|| = ||({diff_vec[0]}, {diff_vec[1]})|| = {length}")

    # 4. State the final answer.
    # By the triangle inequality, ||u-v|| <= ||u|| + ||v|| = 1 + 1 = 2.
    # Since we found an example with length 2, this is the maximum possible.
    print("\nThe largest possible length of such a segment is 2.")


solve()
