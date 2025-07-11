def compare_endomorphisms(rank_A, dim_T):
    """
    Demonstrates the relationship between the number of endomorphisms of a
    semi-abelian variety G and its underlying abelian variety A.

    Args:
      rank_A (int): The rank of the endomorphism ring of the abelian variety A.
      dim_T (int): The dimension of the torus T.
    """

    print(f"Let's assume the rank of the endomorphism ring of A is rank(End(A)) = {rank_A}.")
    print(f"Let's assume the dimension of the toric part T is r = {dim_T}.")
    print("-" * 30)
    
    # The rank of the endomorphism ring of a torus of dimension r is r^2.
    rank_T = dim_T**2
    print(f"The rank of the endomorphism ring of T is rank(End(T)) = r^2 = {dim_T}^2 = {rank_T}.")

    # The rank of End(G) is the sum of the ranks of End(A) and End(T).
    rank_G = rank_A + rank_T

    print("The formula relating the ranks is: rank(End(G)) = rank(End(A)) + rank(End(T)).")
    print("Using the given numbers, we calculate the rank of End(G):")
    # The final equation with each number printed out
    print(f"{rank_G} = {rank_A} + {rank_T}")

    print("-" * 30)
    print(f"Comparing the ranks: rank(End(G)) = {rank_G} and rank(End(A)) = {rank_A}.")

    if rank_G > rank_A:
        print("Conclusion: G has more endomorphisms than A.")
    else: # This happens only when dim_T = 0
        print("Conclusion: G and A have the same number of endomorphisms (since G is A).")
    print("In all cases, rank(End(G)) >= rank(End(A)).")


# --- Demonstration with an example ---
# Example: A is an elliptic curve without complex multiplication (rank=1)
# and T is a 2-dimensional torus (r=2).
example_rank_A = 1
example_dim_T = 2
compare_endomorphisms(example_rank_A, example_dim_T)