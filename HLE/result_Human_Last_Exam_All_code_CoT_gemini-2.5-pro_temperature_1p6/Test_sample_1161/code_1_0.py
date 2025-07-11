import numpy as np

def demonstrate_fortress_problem_failure():
    """
    Demonstrates that a finite number of guards is insufficient to observe
    the exterior of a unit sphere. We use k=4 guards at the vertices of an
    inscribed regular tetrahedron as a concrete example.
    """
    print("Demonstrating the fortress problem for a 3D unit ball.")
    print("We test the case with k=4 guards placed at the vertices of a regular tetrahedron.")
    print("-" * 70)

    # 1. Define 4 guard positions on the unit sphere
    # These are the vertices of a regular tetrahedron.
    c = 1.0 / np.sqrt(3)
    guards = np.array([
        [c, c, c],
        [c, -c, -c],
        [-c, c, -c],
        [-c, -c, c]
    ])

    # 2. Find the least-covered direction
    # By symmetry, the direction least covered by G1, G2, G3 is the
    # direction opposite to G0. Let's verify. We are looking for a direction `u`
    # on the unit sphere that minimizes `max(u . Gi)`. The center of a face
    # is a good candidate for this `u`. The direction to the center of the face
    # defined by G1, G2, G3 is opposite to G0.
    u0 = -guards[0]
    print(f"The chosen least-covered direction u0 is: {np.round(u0, 3)}")


    # 3. Calculate C_min, the maximum coverage in this direction
    # This value must be less than 1 for our argument to hold.
    dot_products = [np.dot(u0, g) for g in guards]
    C_min = np.max(dot_products)

    print(f"The 'coverage' C_min in this direction is max(u0 . Gi) = {C_min:.4f}")
    print("Since C_min < 1, there's a visibility gap in this direction.")

    # 4. Find the range of unseen points
    # Any point X = r*u0 with 1 < r < 1/C_min will be unseen.
    r_limit = 1.0 / C_min
    print(f"This means points X = r*u0 are unseen for 1 < r < {r_limit:.4f}")

    # 5. Show a concrete example of an unseen point
    # Pick a radius r within the unseen range.
    r_example = (1.0 + r_limit) / 2.0
    print(f"\nLet's test a specific point X at radius r = {r_example:.4f}")

    # The condition to be seen is r * C_min >= 1
    unseen_value = r_example * C_min

    print("\nFor X to be seen, the value of 'r * C_min' must be >= 1.")
    print("The final equation for our example point is:")
    # Using f-string to format the output of the "equation"
    print(f"{r_example:.4f} * {C_min:.4f} = {unseen_value:.4f}")

    print(f"\nSince {unseen_value:.4f} is less than 1, this point X is outside the ball but is not seen by any guard.")
    print("This argument holds for any finite number of guards.")


if __name__ == '__main__':
    demonstrate_fortress_problem_failure()