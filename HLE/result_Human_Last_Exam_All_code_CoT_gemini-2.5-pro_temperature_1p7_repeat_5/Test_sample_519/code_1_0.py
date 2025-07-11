def solve_cfg_properties():
    """
    Calculates the properties of the three given categories fibered in groupoids
    and formats the output as requested.
    """

    # --- Properties for X1 ---
    # X1 is the Hilbert scheme of 11 subschemes in A^3.
    # Type: Scheme (S)
    # Properties: separated (s)
    # Dimension: 11 points in a 3-dimensional space.
    dim1 = 3 * 11
    props1 = f"[S, s, {dim1}]"

    # --- Properties for X2 ---
    # X2 is the quotient stack [(A^4 \ V(xy-zw))/C*].
    # Type: Deligne-Mumford stack (DM)
    # Properties: separated (s), irreducible (irr)
    # Dimension: dim(A^4) - dim(C*)
    dim2 = 4 - 1
    props2 = f"[DM, s, irr, {dim2}]"

    # --- Properties for X3 ---
    # X3 is the Picard stack of a genus 7 curve.
    # Type: Algebraic stack (A)
    # Properties: separated (s)
    # Dimension: genus + dim(stabilizer) = g + dim(G_m)
    genus = 7
    dim_stabilizer = 1
    dim3 = genus + dim_stabilizer
    props3 = f"[A, s, {dim3}]"

    # --- Final formatted string ---
    final_answer = f"{props1} {props2} {props3}"
    print(final_answer)

solve_cfg_properties()