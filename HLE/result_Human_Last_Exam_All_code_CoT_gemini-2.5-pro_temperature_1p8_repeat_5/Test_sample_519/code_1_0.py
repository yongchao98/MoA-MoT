import sys

def solve():
    """
    This function analyzes the properties of three categories fibered in groupoids (stacks)
    and prints the results in the specified format.
    """

    # --- Step 1: Analyze X_1 ---
    # X_1 is the Hilbert scheme of 11 points in affine 3-space, Hilb_11(A^3).
    # - Type: Hilbert schemes are representable by schemes. (S)
    # - Separatedness: Hilbert schemes of affine space A^n for n>=2 are generally not separated,
    #   and therefore not universally closed. (no 's', no 'uc')
    # - Irreducibility: Hilb_d(A^n) is known to be reducible for n>=3. Here n=3, d=11. (no 'irr')
    # - Dimension: For Hilb_d(A^n), the dimension is n * d.
    n_1 = 3
    d_1 = 11
    dim_1 = n_1 * d_1
    print(f"Calculating dimension for X1: {n_1} * {d_1} = {dim_1}")
    profile_1 = f"[S,{dim_1}]"

    # --- Step 2: Analyze X_2 ---
    # X_2 = [(A^4 \ V(xy-zw))/C*] with weights (1,4,2,3).
    # - Type: The action of C* on U = A^4 \ V(xy-zw) is free. The GIT quotient of a
    #   quasi-affine variety by a free action of a reductive group like C* is a scheme. (S)
    # - Separatedness: As a good GIT quotient, it is separated. (s)
    # - Universally Closed: The quotient is an open subscheme of the proper weighted projective
    #   space P(1,4,2,3). Since it's a proper open subset, it is not itself proper (uc).
    # - Irreducibility: U is the complement of an irreducible hypersurface in an irreducible
    #   space (A^4), so U is irreducible. The quotient is also irreducible. (irr)
    # - Dimension: The dimension of the quotient is dim(total space) - dim(group).
    dim_U_2 = 4
    dim_G_2 = 1
    dim_2 = dim_U_2 - dim_G_2
    print(f"Calculating dimension for X2: {dim_U_2} - {dim_G_2} = {dim_2}")
    profile_2 = f"[S,s,irr,{dim_2}]"

    # --- Step 3: Analyze X_3 ---
    # X_3 is the Picard stack of a genus 7 curve C0, denoted Pic_{C0}.
    # - Type: The stabilizer of any line bundle is Aut(L) which is isomorphic to C*.
    #   Since the stabilizer group is not finite, it is an Algebraic Stack, not Deligne-Mumford. (A)
    # - Separatedness: It is a G_m-gerbe over the Picard scheme Pic(C0), which is
    #   separated. A gerbe over a separated base is separated. (s)
    # - Universally Closed: The coarse moduli space Pic(C0) is an infinite disjoint union of
    #   proper varieties (the Pic^d(C0)), so it is not quasi-compact, hence not proper (uc).
    # - Irreducibility: Since the coarse space has infinitely many connected components, it is not
    #   irreducible.
    # - Dimension: The dimension of the stack is dim(coarse space) + dim(generic stabilizer).
    genus_3 = 7
    dim_stab_3 = 1
    dim_3 = genus_3 + dim_stab_3
    print(f"Calculating dimension for X3: {genus_3} + {dim_stab_3} = {dim_3}")
    profile_3 = f"[A,s,{dim_3}]"
    
    # --- Final Assembly ---
    # Combine the profiles into the final answer format.
    final_answer = f"{profile_1} {profile_2} {profile_3}"
    print(f"<<<{final_answer}>>>")

solve()