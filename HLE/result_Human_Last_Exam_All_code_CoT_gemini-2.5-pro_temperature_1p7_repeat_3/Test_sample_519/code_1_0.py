def solve_cfg_properties():
    """
    Analyzes three categories fibered in groupoids (CFGs) and determines their properties.
    The properties are:
    - Type: S (scheme), DM (Deligne-Mumford stack), A (Algebraic stack)
    - s: separated
    - uc: universally closed
    - irr: irreducible
    - dim: dimension over C
    """

    # --- Case 1: X_1 ---
    # X_1(S) = subschemes Z in A^3 x S, flat over S with degree 11.
    # This is the Hilbert scheme Hilb^11(A^3).

    # Type: The Hilbert scheme is a scheme.
    x1_type = 'S'
    # Separatedness: Hilbert schemes of subschemes of quasi-projective schemes are separated.
    x1_s = 's'
    # Universally Closed: Hilb^n(A^m) is not proper (not compact) for m > 0, as points can go to infinity.
    x1_uc = False
    # Irreducibility: Hilb^n(A^m) is reducible for m >= 3 and n >= 4. Here n=11, m=3.
    x1_irr = False
    # Dimension: The dimension of the component of n distinct points in A^m is n*m.
    # This component is known to have the maximal dimension.
    n_1 = 11
    m_1 = 3
    x1_dim = n_1 * m_1
    
    # Construct profile string for X1
    profile1_list = [x1_type, x1_s]
    if x1_uc: profile1_list.append('uc')
    if x1_irr: profile1_list.append('irr')
    profile1_list.append(str(x1_dim))
    profile1_str = f"[{','.join(profile1_list)}]"


    # --- Case 2: X_2 ---
    # X_2 = [(A^4 \ V(xy-zw))/C*] with weights (1,4,2,3).
    # This is a quotient stack.

    # Type: The stabilizer of any point is trivial. For a point (x,y,z,w), if x is nonzero,
    # the weight 1 is active, so gcd of weights is 1. If x is zero, then zw cannot be zero,
    # so weights 2 and 3 are active, and gcd(2,3)=1. Trivial (finite) stabilizers mean it's a DM stack.
    x2_type = 'DM'
    # Separatedness: A quotient stack of a separated scheme (A^4 \ V(f)) by an affine group (C*) is separated.
    x2_s = 's'
    # Universally Closed: The space A^4 \ V(xy-zw) is not proper, so the quotient stack is not proper.
    x2_uc = False
    # Irreducibility: The space U = A^4 \ V(xy-zw) is irreducible because xy-zw is an irreducible polynomial.
    # The quotient of an irreducible space by a connected group (C*) is irreducible.
    x2_irr = 'irr'
    # Dimension: The dimension of a quotient stack [U/G] is dim(U) - dim(G).
    dim_U_2 = 4  # Dimension of A^4 \ V(f) is 4
    dim_G_2 = 1  # Dimension of C*
    x2_dim = dim_U_2 - dim_G_2
    
    # Construct profile string for X2
    profile2_list = [x2_type, x2_s]
    if x2_uc: profile2_list.append('uc')
    if x2_irr: profile2_list.append(x2_irr)
    profile2_list.append(str(x2_dim))
    profile2_str = f"[{','.join(profile2_list)}]"

    # --- Case 3: X_3 ---
    # X_3(S) = line bundles L on S x C_0, for a genus 7 curve C_0.
    # This is the Picard stack Pic(C_0).

    # Type: The automorphism group of any line bundle on a curve is C*, which is not finite.
    # So this is an Algebraic stack, not Deligne-Mumford.
    x3_type = 'A'
    # Separatedness: The Picard stack is a standard example of a separated algebraic stack.
    x3_s = 's'
    # Universally Closed: It consists of an infinite disjoint union of components (one for each degree d in Z).
    # A scheme or stack with infinitely many components is not of finite type, hence not proper.
    x3_uc = False
    # Irreducibility: Since it has infinitely many connected components, it is not irreducible.
    x3_irr = False
    # Dimension: The dimension of a stack is the max dimension of its components.
    # Each component Pic^d(C_0) is a B(C*)-gerbe over the Jacobian J(C_0), so its dimension is dim(J(C_0)) - dim(C*).
    # dim(J(C_0)) = genus.
    g_3 = 7
    dim_G_3 = 1  # Dimension of C*
    x3_dim = g_3 - dim_G_3

    # Construct profile string for X3
    profile3_list = [x3_type, x3_s]
    if x3_uc: profile3_list.append('uc')
    if x3_irr: profile3_list.append('irr')
    profile3_list.append(str(x3_dim))
    profile3_str = f"[{','.join(profile3_list)}]"

    # Combine and print the final result
    final_answer = f"{profile1_str} {profile2_str} {profile3_str}"
    print(final_answer)

solve_cfg_properties()