def solve():
    """
    This function calculates the properties of the three given stacks and prints them
    in the required format.
    """

    # --- Analysis of X1 ---
    # X1 is the Hilbert scheme of subschemes of degree 11 in A^3.
    # Type: Scheme (S)
    # Properties: separated (s)
    # Dimension: For the main component of 0-dimensional subschemes of length d in A^n,
    # the dimension is n * d.
    n1 = 3
    d1 = 11
    dim1 = n1 * d1
    x1_profile = f"[S,s,{dim1}]"

    # --- Analysis of X2 ---
    # X2 is the quotient stack [ (A^4 \ V(xy-zw)) / C* ].
    # Type: Deligne-Mumford stack (DM) due to finite stabilizers.
    # Properties: separated (s), irreducible (irr).
    # Dimension: dim(total space) - dim(group).
    dim_space = 4
    dim_group = 1
    dim2 = dim_space - dim_group
    x2_profile = f"[DM,s,irr,{dim2}]"

    # --- Analysis of X3 ---
    # X3 is the Picard stack of a genus 7 curve.
    # Type: Algebraic stack (A) due to C* stabilizers.
    # Properties: separated (s).
    # Dimension: Equal to the genus of the curve.
    g = 7
    dim3 = g
    x3_profile = f"[A,s,{dim3}]"

    # --- Final Answer ---
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve()