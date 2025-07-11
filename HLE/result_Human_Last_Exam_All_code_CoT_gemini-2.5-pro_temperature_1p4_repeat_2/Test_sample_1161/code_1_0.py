def solve_fortress_problem_for_sphere():
    """
    Calculates the minimum number of guards needed to observe the exterior of a unit ball in 3D space,
    with guards placed on the ball's surface.
    """

    # The problem is a specific instance of the "fortress problem" or "exterior illumination problem".
    # For a convex body in a d-dimensional space, there is a known result for the minimum
    # number of guards required.
    # A theorem by R. Blind and G. Blind states that for a Euclidean ball in R^d,
    # the minimum number of guards is 2d.

    # In our case, the space is R^3, so the dimension d is 3.
    dimension = 3

    # The formula for the minimum number of guards is 2 * d.
    num_guards = 2 * dimension

    print("The fortress problem for a unit ball in 3D space asks for the minimum number of guards on the surface to see the entire exterior.")
    print("According to a theorem by Blind & Blind, the number of guards for a ball in d-dimensional space is 2d.")
    print(f"For a 3D space, the dimension d is {dimension}.")
    print(f"The final calculation is: 2 * {dimension} = {num_guards}")
    print(f"Therefore, the minimum number of guards necessary is {num_guards}.")

solve_fortress_problem_for_sphere()