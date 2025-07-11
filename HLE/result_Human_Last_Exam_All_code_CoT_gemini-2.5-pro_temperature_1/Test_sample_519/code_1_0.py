def solve_cfgs():
    """
    This function analyzes three categories fibered in groupoids (CFGs)
    and determines their properties, printing the result in the specified format.
    """

    # Analysis for X1: Hilbert scheme of 11 points in A^3
    # Type: Scheme (S)
    # Properties: separated (s). Not universally closed. Not irreducible.
    # Dimension: n * d for the principal component.
    n_x1 = 3  # for A^3
    d_x1 = 11 # for degree 11
    dim_x1 = n_x1 * d_x1
    profile_x1 = f"[S,s,{dim_x1}]"

    # Analysis for X2: Quotient stack [ (A^4 \ V(xy-zw)) / C* ]
    # Type: Deligne-Mumford stack (DM), as the action is free.
    # Properties: separated (s), irreducible (irr). Not universally closed.
    # Dimension: dim(space) - dim(group).
    dim_space_x2 = 4 # dim(A^4)
    dim_group_x2 = 1 # dim(C*)
    dim_x2 = dim_space_x2 - dim_group_x2
    profile_x2 = f"[DM,s,irr,{dim_x2}]"

    # Analysis for X3: Picard stack of a genus 7 curve
    # Type: Algebraic stack (A), as stabilizers are C*.
    # Properties: separated (s). Not universally closed. Not irreducible.
    # Dimension: genus of the curve.
    dim_x3 = 7 # given genus
    profile_x3 = f"[A,s,{dim_x3}]"

    # Combine the profiles and print the final answer
    final_answer = f"{profile_x1} {profile_x2} {profile_x3}"
    print(final_answer)

solve_cfgs()