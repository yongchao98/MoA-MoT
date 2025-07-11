def solve():
    """
    Analyzes the seven statements and provides the correctness of each.
    a) The quotient topology on Omega_g^rt / GL_g(Z) is Hausdorff. (False)
    b) The barycentric subdivision of the perfect cone decomposition of Omega_4^rt consists only of simplicial cones. (True)
    c) For any compact set K in Omega_g^rt there exists a set of finitely many cones in the second Voronoi compactification such that their union contains K. (True)
    d) The number of orbits of maximal cones in the perfect cone decomposition of Omega_7^rt is 33. (True)
    e) The GL_g(Z)-stabilizer of a cone sigma in the central cone decomposition is finite if sigma intersects Omega_g. (True)
    f) The number of orbits of cones in the Second Voronoi decomposition of Omega_5^rt is 222. (True)
    g) Given cones tau, sigma in the perfect cone decomposition with tau a face of sigma, one has that their corresponding groups of GL_g(Z)-stabilizers satisfy Stab(tau) is a subset of Stab(sigma). (False)
    """
    
    # The reasoning for each statement is provided in the text above.
    # a) N - Not Hausdorff at the origin.
    # b) Y - Barycentric subdivision is always simplicial.
    # c) Y - This is the property of local finiteness.
    # d) Y - Known result by Jaquet-Chiffelle.
    # e) Y - Stabilizer of a cone of the first kind is finite.
    # f) Y - Known result by Dutour, Schürmann, Vallentin.
    # g) N - An element stabilizing a face can permute the cones containing it.
    
    answer = "NYYYYTN"
    print(answer)

solve()