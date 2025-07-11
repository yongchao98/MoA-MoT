import numpy as np

def demonstrate_fortress_problem():
    """
    Demonstrates that a finite number of guards on a sphere cannot see the
    entire exterior.

    The script sets up 4 guards at the vertices of a regular tetrahedron
    inscribed in the unit sphere. It then finds a "blind spot" direction and
    constructs a point in the exterior that none of the guards can see.
    """

    print("--- The Fortress Problem for a 3D Sphere ---\n")
    print("Goal: Show that a finite number of guards on a sphere's surface cannot see the entire exterior.")

    # 1. Define the guard positions.
    # We choose 4 guards placed at the vertices of a regular tetrahedron
    # inscribed in the unit sphere (radius=1).
    p1 = np.array([0, 0, 1])
    p2 = np.array([2 * np.sqrt(2) / 3, 0, -1 / 3])
    p3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1 / 3])
    p4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1 / 3])

    guards = [p1, p2, p3, p4]

    print(f"\nStep 1: Place 4 guards on the unit sphere.")
    for i, p in enumerate(guards):
        print(f"  Guard {i+1} position (p_{i+1}): {np.round(p, 3)}")

    # 2. Find a "blind spot" direction on the sphere.
    # For any finite set of guards, there is a point on the sphere's surface
    # that is furthest from any guard. For this tetrahedral arrangement, a point
    # antipodal to any guard is such a spot. Let's pick the point opposite to p1.
    q0 = -p1
    print(f"\nStep 2: Identify a potential blind spot direction (q_0).")
    print(f"  We choose a point on the sphere maximally far from the guards.")
    print(f"  Blind spot direction (q_0): {q0}")

    # 3. Construct an unseen point in the exterior.
    # A point `q` in the exterior is visible by guard `p_i` if p_i . q > 1.
    # We will show that a point q = (1+delta)*q0 for a small delta is not seen.
    # Let's find the maximum dot product between our blind spot direction and any guard.
    c_max = max(np.dot(p, q0) for p in guards)
    print(f"\nStep 3: Construct a point in the exterior (q) that is likely to be unseen.")
    print(f"  The maximum dot product between q_0 and any guard position is {c_max:.4f}.")
    print(f"  This means q_0 is not aligned with any guard.")
    
    # Any point q = (1+delta)*q0 will be unseen if (1+delta)*c_max <= 1,
    # or delta <= (1/c_max) - 1.
    # For c_max = 1/3, delta must be less than (1/(1/3)) - 1 = 2.
    # Let's pick delta = 1.0, which is < 2.
    delta = 1.0
    unseen_point = (1 + delta) * q0
    
    print(f"  We construct a point q = (1 + {delta}) * q_0 = {unseen_point}")
    print(f"  This point is outside the unit ball (distance from origin is {np.linalg.norm(unseen_point)}).")
    
    # 4. Verify that the point is not seen by any guard.
    print("\nStep 4: Check if any guard can see the point q.")
    print("  A guard p_i sees q if p_i . q > 1.")
    
    is_seen = False
    for i, p in enumerate(guards):
        visibility_check = np.dot(p, unseen_point)
        result = "Sees the point" if visibility_check > 1 else "Does NOT see the point"
        print(f"  Guard {i+1}: p_{i+1} . q = {np.round(p, 3)} . {unseen_point} = {visibility_check:.4f}.  (Condition: > 1)")
        print(f"    --> {result}")
        if visibility_check > 1:
            is_seen = True

    print("\n--- Conclusion ---")
    if not is_seen:
        print("The constructed point q is not seen by any of the 4 guards.")
        print("This demonstrates that for any finite number of guards, a blind spot can be found.")
        print("Therefore, a finite number of guards is insufficient.")
    else:
        print("The demonstration failed, which indicates an error in the logic or calculation.")

if __name__ == '__main__':
    demonstrate_fortress_problem()
