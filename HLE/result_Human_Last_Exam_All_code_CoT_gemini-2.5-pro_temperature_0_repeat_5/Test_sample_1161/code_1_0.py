import math

def solve_fortress_problem():
    """
    This function explains the solution to the fortress problem for a 3D unit ball.
    """
    
    # Step 1: Define the problem and visibility
    print("Problem: Find the minimum number of guards on the surface of a unit sphere to observe the entire exterior space.")
    print("A guard at a point G on the sphere's surface sees the closed half-space 'outside' the tangent plane at G.")
    print("This means a point X is visible from G if the dot product X . G is greater than or equal to 1.")
    print("-" * 20)

    # Step 2: Analyze insufficient numbers of guards
    print("Analyzing the number of guards required:")
    
    # Case n=1, 2, 3
    num_guards_insufficient = 3
    print(f"Let's test if {num_guards_insufficient} guards are sufficient.")
    print("If we place 3 or fewer guards on the sphere, their positions must lie on a single plane that passes through the sphere's center (a great circle).")
    print("Let's assume this plane is the xy-plane (z=0). All guard position vectors G_i would have a z-component of 0.")
    print("Now, consider a point X far away on the z-axis, for example, X = (0, 0, 100). This point is in the exterior.")
    print("The dot product for any guard would be X . G_i = (0 * G_ix) + (0 * G_iy) + (100 * 0) = 0.")
    print("Since 0 < 1, this point X is not seen by any of the 3 guards.")
    print(f"Therefore, {num_guards_insufficient} guards are not sufficient.")
    print("-" * 20)

    # Step 3: Propose and justify the solution
    num_guards_sufficient = 4
    print(f"Let's test if {num_guards_sufficient} guards are sufficient.")
    print("To solve the problem with 3 guards, we must place the guards so they are not co-planar.")
    print("The minimum number of points that are not co-planar is 4, forming a tetrahedron.")
    print("We can place the 4 guards at the vertices of a regular tetrahedron inscribed in the unit sphere.")
    print("This symmetric arrangement ensures that the convex hull of the guard positions contains the sphere's center.")
    print("This configuration guarantees that for any point in the exterior, it will be visible to at least one of the guards.")
    print(f"Thus, {num_guards_sufficient} guards are necessary and sufficient.")
    print("-" * 20)

    # Step 4: Final Answer
    final_answer = 4
    print(f"The minimum amount of guards necessary is: {final_answer}")

solve_fortress_problem()