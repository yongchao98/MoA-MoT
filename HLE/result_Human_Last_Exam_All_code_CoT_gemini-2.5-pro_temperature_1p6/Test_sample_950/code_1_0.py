import sys

def solve_task():
    """
    This script calculates the rank of the torsion subgroup of the integral
    cohomology ring of the space of 3-subspaces of R^5.
    """
    
    # Step 1: Identify the mathematical space.
    # The space of k-subspaces in R^n is the real Grassmannian manifold Gr(k, n).
    # For this problem, we have k=3 and n=5.
    k = 3
    n = 5
    
    print(f"The space in question is the real Grassmannian manifold Gr({k}, {n}).")

    # Step 2: Use a known property of Grassmannians.
    # There is a diffeomorphism Gr(k, n) is congruent to Gr(n-k, n).
    # So, Gr(3, 5) is diffeomorphic to Gr(5-3, 5) = Gr(2, 5).
    # We will analyze the cohomology of Gr(2, 5) as it has the same structure.
    print(f"This space is diffeomorphic to Gr({n-k}, {n}), so we study the cohomology of Gr({n-k}, {n}).")

    # Step 3: State the known structure of the integral cohomology ring.
    # The integral cohomology ring, H*(Gr(k, n); Z), is the direct sum of
    # cohomology groups H^i(Gr(k, n); Z). Unlike complex Grassmannians, the
    # integral cohomology of real Grassmannians can have 2-torsion.
    # The results for Gr(2, 5) are known from the literature in algebraic topology.
    print("\nThe integral cohomology groups H^i(Gr(2, 5); Z) with non-trivial torsion are:")
    print("H^3(Gr(2, 5); Z) is isomorphic to Z_2")
    print("H^4(Gr(2, 5); Z) contains a torsion part of Z_2")
    print("H^5(Gr(2, 5); Z) is isomorphic to Z_2")
    
    # Step 4: Identify the torsion subgroup of the cohomology ring.
    # The torsion subgroup of the entire ring is the direct sum of the torsion
    # parts of each H^i. Let T be this subgroup.
    # T is isomorphic to Z_2 + Z_2 + Z_2.
    print("\nThe torsion subgroup T of the whole cohomology ring H*(Gr(2, 5); Z) is the direct sum of these torsion parts.")
    print("Therefore, T is isomorphic to Z_2 + Z_2 + Z_2, or (Z_2)^3.")

    # Step 5: Compute the rank of this torsion subgroup.
    # The 'rank' of a finitely generated abelian group is the number of copies of Z
    # in its structure decomposition (i.e., the rank of its free part).
    # A torsion group consists entirely of elements of finite order and thus has no free part (or a free part of rank 0).
    final_rank = 0
    print("\nThe rank of an abelian group is the rank of its free component.")
    print("Since T is a torsion group, its free component is trivial.")
    
    print("\nFinal calculation:")
    print(f"rank(Torsion(H*(Gr({k}, {n}); Z))) = {final_rank}")

solve_task()