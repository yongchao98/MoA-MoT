def identify_manifold():
    """
    Identifies the three-manifold from the given Heegaard diagram and prints its properties.
    """
    p, q, r = 2, 3, 4

    print("The Heegaard diagram represents the Brieskorn sphere Sigma(p, q, r).")
    print(f"For this diagram, the parameters are (p, q, r) = ({p}, {q}, {r}).")
    print("So the manifold is Sigma(2, 3, 4).\n")

    print("This is a spherical 3-manifold, also known as a prism manifold.")
    print("Its fundamental group, pi_1(Sigma(2, 3, 4)), is the binary tetrahedral group of order 24.\n")

    print("A standard presentation for the fundamental group of Sigma(p, q, r) is:")
    print("<a, b, c | a^p = b^q, b^q = c^r, abc = 1>\n")

    print("For this specific manifold, the defining relations are:")
    # The prompt requires printing each number in the final equation.
    a_exponent = p
    b_exponent = q
    c_exponent = r
    
    print(f"a^{a_exponent} = b^{b_exponent}")
    print(f"b^{b_exponent} = c^{c_exponent}")
    print("a * b * c = 1")

identify_manifold()