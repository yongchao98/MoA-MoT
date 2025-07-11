def solve():
    """
    Determines the correctness of seven statements related to the rational closure of the cone of positive definite matrices.

    a) The quotient topology on Omega_g^rt / GL_g(Z) is Hausdorff. (True)
    b) The barycentric subdivision of the perfect cone decomposition of Omega_4^rt consists only of simplicial cones. (True)
    c) For any compact set K in Omega_g^rt there exists a set of finitely many cones in the second Voronoi compactification such that their union contains K. (True)
    d) The number of orbits of maximal cones in the perfect cone decomposition of Omega_7^rt is 33. (True)
    e) The GL_g(Z)-stabilizer of a cone sigma in the central cone decomposition is finite if sigma intersects Omega_g. (True)
    f) The number of orbits of cones in the Second Voronoi decomposition of Omega_5^rt is 222. (False, 222 is for g=7, for g=5 it's 20)
    g) Given cones tau, sigma in the perfect cone decomposition with tau a face of sigma, one has that their corresponding groups of GL_g(Z)-stabilizers satisfy Stab(tau) is a subset of Stab(sigma). (False, the inclusion is reversed)
    """
    
    # Based on the analysis, the answers are:
    # a: Y
    # b: Y
    # c: Y
    # d: Y
    # e: Y
    # f: N
    # g: N
    
    answer_string = "YYYYYNN"
    print(answer_string)

solve()