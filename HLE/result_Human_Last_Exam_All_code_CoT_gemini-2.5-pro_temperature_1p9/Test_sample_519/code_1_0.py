def solve_cfgs():
    """
    Analyzes three categories fibered in groupoids and determines their properties.
    """

    # Properties for X_1(S) = subschemes Z in A^3 x S, flat over S with degree 11
    # This is the Hilbert Scheme Hilb^11(A^3)
    type_1 = "S"         # It's a Scheme
    sep_1 = "s"          # It's separated
    # Not universally closed, not irreducible
    dim_1 = 3 * 11       # Dimension is n*d = 3*11 = 33
    profile_1 = f"[{type_1},{sep_1},{dim_1}]"

    # Properties for X_2 = [ (A^4 \ V(xy-zw)) / C* ] with weights (1,4,2,3)
    type_2 = "DM"        # It's a Deligne-Mumford stack (finite stabilizers)
    sep_2 = "s"          # It's separated
    irr_2 = "irr"        # It's irreducible
    # Not universally closed
    dim_2 = 4 - 1        # Dimension is dim(space) - dim(group) = 4 - 1 = 3
    profile_2 = f"[{type_2},{sep_2},{irr_2},{dim_2}]"

    # Properties for X_3(S) = line bundles L on S x C_0, where C_0 is a genus 7 curve
    # This is the Picard Scheme Pic(C_0)
    type_3 = "S"         # It's a Scheme
    sep_3 = "s"          # It's separated
    # Not universally closed, not irreducible
    dim_3 = 7            # Dimension is the genus of the curve, g=7
    profile_3 = f"[{type_3},{sep_3},{dim_3}]"

    # Combine the profiles into a single string
    final_answer = f"{profile_1} {profile_2} {profile_3}"

    print(final_answer)

solve_cfgs()