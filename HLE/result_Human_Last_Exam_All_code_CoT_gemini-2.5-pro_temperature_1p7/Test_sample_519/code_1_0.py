def solve():
    """
    This function computes the properties of the three given CFGs and formats the output.
    """

    # --- Properties for X1 ---
    # X1 is the Hilbert scheme Hilb_11(A^3).
    # Type: Scheme (S)
    # Separated: Yes (s)
    # Universally Closed: No
    # Irreducible: No
    # Dimension: n * d = 3 * 11
    dim1_n = 3
    dim1_d = 11
    dim1 = dim1_n * dim1_d
    profile1 = f"[S,s,{dim1}]"

    # --- Properties for X2 ---
    # X2 is the quotient stack [(A^4 \ V(xy-zw))/C*].
    # Type: Scheme (S)
    # Separated: Yes (s)
    # Universally Closed: No
    # Irreducible: Yes (irr)
    # Dimension: dim(A^4) - dim(C*) = 4 - 1
    dim2_space = 4
    dim2_group = 1
    dim2 = dim2_space - dim2_group
    profile2 = f"[S,s,irr,{dim2}]"

    # --- Properties for X3 ---
    # X3 is the Picard stack Pic(C_0) for a genus 7 curve C_0.
    # Type: Algebraic Stack (A)
    # Separated: Yes (s)
    # Universally Closed: No
    # Irreducible: No
    # Dimension: genus g = 7
    dim3_g = 7
    dim3 = dim3_g
    profile3 = f"[A,s,{dim3}]"

    # --- Combine and print the results ---
    final_answer = f"{profile1} {profile2} {profile3}"
    print(final_answer)

solve()