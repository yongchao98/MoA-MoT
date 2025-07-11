def solve_statements():
    """
    Evaluates the correctness of the seven statements.

    The final answer is formatted as a single string of 'Y' (Yes/True) and 'N' (No/False)
    characters, corresponding to statements a) through g).

    The numbers mentioned in the statements are:
    b) g = 4
    d) g = 7, orbits = 33
    f) g = 5, orbits = 222
    """

    # a) The quotient topology on Omega_g^rt / GL_g(Z) is Hausdorff. -> Y
    # b) The barycentric subdivision of the perfect cone decomposition of Omega_4^rt consists only of simplicial cones. -> Y
    # c) For any compact set K in Omega_g^rt, it is contained in a finite union of cones. -> N
    # d) The number of orbits of maximal cones for g=7 is 33. -> Y
    # e) The GL_g(Z)-stabilizer of a cone sigma is finite if sigma intersects Omega_g. -> N
    # f) The number of orbits of cones for g=5 is 222. -> N
    # g) For a face tau of a cone sigma, Stab(tau) is a subset of Stab(sigma). -> Y

    final_answer = "YYNYNNY"
    print(final_answer)

solve_statements()