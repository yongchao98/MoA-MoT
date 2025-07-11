def solve_cfgs():
    """
    Determines and prints the properties of three specified categories fibered in groupoids.

    The properties for each category are determined based on standard results in algebraic geometry and stack theory.
    - X_1 is the Hilbert scheme of 11 points in A^3.
    - X_2 is a quotient stack [ (A^4 \ V(xy-zw)) / C* ].
    - X_3 is the Picard stack of a genus 7 curve.
    """

    # Profile for X_1: Hilb^11(A^3)
    # [Type=Scheme, Separated=yes, Univ.Closed=no, Irreducible=no, Dimension=3*11=33]
    profile_x1 = "[S,s,33]"

    # Profile for X_2: [ (A^4 \ V(xy-zw)) / C* ] with weights (1,4,2,3)
    # [Type=DM stack, Separated=yes, Univ.Closed=no, Irreducible=yes, Dimension=4-1=3]
    profile_x2 = "[DM,s,irr,3]"

    # Profile for X_3: Picard stack of a genus 7 curve
    # [Type=Algebraic stack, Separated=no, Univ.Closed=no, Irreducible=no, Dimension=g=7]
    profile_x3 = "[A,7]"

    final_answer = f"{profile_x1} {profile_x2} {profile_x3}"
    print(final_answer)

solve_cfgs()