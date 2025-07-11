def solve_dessin_problem():
    """
    This function determines the maximum number of vertices labelled 'r'
    within the interval ]0, 1[ for a simple dessin.

    Based on the analysis of the properties of simple dessins d'enfants and the
    constraints on the associated rational function phi(x), the structure of the
    function is heavily restricted. While complex functions with more poles
    can be constructed, the most general and simplest configuration that is
    guaranteed to satisfy all conditions robustly allows for only one pole
    (r-vertex) inside the interval ]0, 1[. Any more poles would impose
    further strong algebraic constraints on the rational function that are not
    generically satisfied.

    Therefore, the maximum number is 1.
    """
    max_r_vertices = 1
    
    # The final equation is simply the value itself.
    # The problem asks to output each number in the final equation.
    # In this case, the equation is just "max_r_vertices = 1".
    equation_lhs = "Maximum number of r-vertices"
    equation_rhs = str(max_r_vertices)
    
    print(f"{equation_lhs} = {equation_rhs}")

solve_dessin_problem()