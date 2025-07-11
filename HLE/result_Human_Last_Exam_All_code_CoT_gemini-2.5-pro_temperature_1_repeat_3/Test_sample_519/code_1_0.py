def solve_cfgs():
    """
    This function provides the properties of three categories fibered in groupoids.
    """
    # Properties for X1: Hilb_A^3^11
    profile_x1 = "[S,s,irr,33]"

    # Properties for X2: [ (A^4 \ V(xy-zw)) / C* ]
    profile_x2 = "[S,s,irr,3]"

    # Properties for X3: Pic(C_0) for a genus 7 curve C_0
    profile_x3 = "[S,s,7]"

    # The final answer is the space-separated list of these profiles.
    final_answer = f"{profile_x1} {profile_x2} {profile_x3}"

    print(final_answer)

solve_cfgs()