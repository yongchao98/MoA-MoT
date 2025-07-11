def solve_cfgs():
    """
    This function analyzes the properties of the three given categories fibered
    in groupoids (CFGs) and prints the result in the specified format.
    """

    # --- Analysis for X_1 ---
    # X_1 is the Hilbert scheme of 11 subschemes in A^3.
    # It is a scheme (S).
    # It is separated (s).
    # It is not universally closed as A^3 is not proper.
    # It is reducible for dimension >= 3.
    # Dimension = degree * ambient_dimension
    x1_degree = 11
    x1_ambient_dim = 3
    x1_dim = x1_degree * x1_ambient_dim
    # Format: [Type, separatedness, dimension]
    x1_profile = f"[S, s, {x1_dim}]"

    # --- Analysis for X_2 ---
    # X_2 is the quotient stack [ (A^4 \ V(xy-zw)) / C* ].
    # The C* action with weights (1,4,2,3) on A^4 \ V(xy-zw) is free,
    # making the quotient a scheme (S).
    # The quotient of a separated scheme by a group with closed orbits is separated (s).
    # The quotient is not universally closed.
    # It is the quotient of an irreducible space, so it is irreducible (irr).
    # Dimension = dim(space) - dim(group)
    x2_space_dim = 4
    x2_group_dim = 1
    x2_dim = x2_space_dim - x2_group_dim
    # Format: [Type, separatedness, irreducibility, dimension]
    x2_profile = f"[S, s, irr, {x2_dim}]"

    # --- Analysis for X_3 ---
    # X_3 is the Picard stack of a genus 7 curve.
    # It has C* stabilizers, so it is an Algebraic Stack (A).
    # It is separated (s).
    # It is not universally closed (it is a Gm-gerbe).
    # It is not irreducible (it's a disjoint union over degrees).
    # Dimension = genus of the curve
    x3_genus = 7
    x3_dim = x3_genus
    # Format: [Type, separatedness, dimension]
    x3_profile = f"[A, s, {x3_dim}]"

    # --- Combine and print the results ---
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve_cfgs()