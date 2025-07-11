def solve_cfg_properties():
    """
    Calculates the properties of the three given categories fibered in groupoids
    and formats the output as specified.
    """

    # Properties for X_1: The Hilbert scheme Hilb^11(A^3)
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally Closed: no
    # Irreducible: no
    # Dimension = n * d = 3 * 11 = 33
    dim_X1 = 3 * 11
    profile_X1 = f"[S, s, {dim_X1}]"

    # Properties for X_2: The quotient stack [ (A^4 \ V(xy-zw)) / C* ]
    # Type: Algebraic Space (S), because the C* action is free.
    # Separated: yes (s), because the action is proper.
    # Universally Closed: no
    # Irreducible: yes (irr)
    # Dimension = dim(space) - dim(group) = 4 - 1 = 3
    dim_X2 = 4 - 1
    profile_X2 = f"[S, s, irr, {dim_X2}]"

    # Properties for X_3: The Picard stack Pic_C0 for a genus 7 curve C0
    # Type: Algebraic Stack (A), because stabilizers are C*.
    # Separated: yes (s)
    # Universally Closed: no
    # Irreducible: no
    # Dimension = genus + dim(stabilizer) = 7 + 1 = 8
    g = 7
    dim_X3 = g + 1
    profile_X3 = f"[A, s, {dim_X3}]"

    # Combine the profiles into a single string
    final_answer = f"{profile_X1} {profile_X2} {profile_X3}"
    
    print(final_answer)

solve_cfg_properties()