def solve_cfgs():
    """
    This function determines the properties for three given categories fibered
    in groupoids and prints them in the specified format.
    """

    # Properties are encoded as follows:
    # S: Scheme or Algebraic Space
    # DM: Deligne-Mumford stack
    # A: Algebraic stack
    # s: separated
    # uc: universally closed
    # irr: irreducible
    # dim: dimension over C

    # --- X_1 ---
    # Hilb_11(A^3)
    # Type: Scheme [S]
    # Irreducible: Yes, by Fogarty's theorem for Hilb_d(A^n), n>=3. [irr]
    # Separated: No. Not 's'.
    # Universally Closed: No, since not separated. Not 'uc'.
    # Dimension: n*d = 3*11=33
    x1_dim = 3 * 11
    x1_profile = f"[S,irr,{x1_dim}]"

    # --- X_2 ---
    # [(A^4 \ V(xy-zw)) / C*] with weights (1,4,2,3)
    # Type: The stack has trivial geometric stabilizers, so it is an algebraic space. [S]
    #       The coarse moduli space is a quasi-projective scheme.
    # Separated: Yes, as a quasi-projective scheme. [s]
    # Irreducible: Yes, as a quotient of an irreducible space. [irr]
    # Universally Closed: No, as it's not projective. Not 'uc'.
    # Dimension: dim(A^4) - dim(C*) = 4 - 1 = 3
    x2_dim = 4 - 1
    x2_profile = f"[S,s,irr,{x2_dim}]"

    # --- X_3 ---
    # Pic(C_0) for a genus 7 curve C_0
    # Type: Scheme [S]
    # Separated: Yes, Picard schemes are separated. [s]
    # Irreducible: No, it is a disjoint union over degrees d in Z.
    # Universally Closed: No, it is not of finite type. Not 'uc'.
    # Dimension: The dimension of each component is the genus, g=7.
    x3_dim = 7
    x3_profile = f"[S,s,{x3_dim}]"

    # Combine and print the results
    final_result = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_result)

solve_cfgs()