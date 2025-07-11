def solve():
    """
    This function calculates and prints the properties of the three given CFGs.
    """

    # Properties for X1: Hilb_11(A^3)
    # Dimension calculation: n * d
    dim_x1 = 3 * 11
    # Format: [Type, separated, dimension]
    # Not universally closed, not irreducible.
    x1_profile = f"[S,s,{dim_x1}]"

    # Properties for X2: [(\A^4 \ V(xy-zw)) / C*]
    # Dimension calculation: dim(space) - dim(group)
    dim_x2 = 4 - 1
    # Format: [Type, separated, irreducible, dimension]
    # Not universally closed.
    x2_profile = f"[DM,s,irr,{dim_x2}]"

    # Properties for X3: Pic(C_0) for a genus 7 curve
    # Dimension calculation: genus + 1
    dim_x3 = 7 + 1
    # Format: [Type, separated, universally closed, dimension]
    # Not irreducible.
    x3_profile = f"[A,s,uc,{dim_x3}]"

    # Combine the profiles into a single string and print
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve()