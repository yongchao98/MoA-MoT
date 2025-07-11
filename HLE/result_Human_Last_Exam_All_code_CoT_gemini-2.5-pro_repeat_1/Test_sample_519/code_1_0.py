def solve_cfgs():
    """
    This function calculates and prints the properties of the three given
    categories fibered in groupoids (CFGs).
    """

    # Properties for X1: Hilb^11(A^3)
    x1_profile = "[S,s,33]"

    # Properties for X2: [ (A^4 \ V(xy-zw)) / C* ]
    x2_profile = "[DM,s,irr,3]"

    # Properties for X3: Pic_{C_0} for a genus 7 curve C_0
    x3_profile = "[A,s,7]"

    # Combine the profiles into a single string as requested
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"

    print(final_answer)

solve_cfgs()