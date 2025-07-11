import math

def solve_grassmannian_torsion_rank():
    """
    This script calculates the rank of the torsion subgroup of the integral 
    cohomology ring of the space of 3-subspaces of R^5.

    This space is the Grassmannian Gr(3, 5).
    """
    n = 5
    k = 3
    n_minus_k = n - k

    # Helper function for SO(m) rank
    def rank_so(m):
        return math.floor(m / 2)

    print("This program calculates the rank of the torsion subgroup of H*(Gr(3, 5); Z).")
    print(f"The space Gr({k}, {n}) is represented as the homogeneous space G/H = SO({n}) / (SO({k}) x SO({n_minus_k})).")
    
    print("\n--- Step 1: Calculate the ranks of the Lie groups G and H ---")
    
    # Rank of G = SO(n)
    rank_g = rank_so(n)
    print(f"The rank of G = SO({n}) is floor({n}/2).")
    print(f"rank(G) = {rank_g}")

    # Rank of H = SO(k) x SO(n-k)
    rank_h_1 = rank_so(k)
    rank_h_2 = rank_so(n_minus_k)
    rank_h = rank_h_1 + rank_h_2
    print(f"The rank of H = SO({k}) x SO({n_minus_k}) is rank(SO({k})) + rank(SO({n_minus_k})).")
    print(f"rank(SO({k})) = floor({k}/2) = {rank_h_1}")
    print(f"rank(SO({n_minus_k})) = floor({n_minus_k}/2) = {rank_h_2}")
    print(f"rank(H) = {rank_h_1} + {rank_h_2} = {rank_h}")
    
    print("\n--- Step 2: Apply Borel's theorem ---")
    print(f"Comparing the ranks: rank(G) = {rank_g} and rank(H) = {rank_h}.")
    if rank_g == rank_h:
        print("Since rank(G) = rank(H), a theorem by A. Borel implies that the integral cohomology ring H*(Gr(3, 5); Z) is torsion-free.")
        print("This means the torsion subgroup is the trivial group {0}.")
    else:
        # This case is not met, but included for completeness
        print("Since rank(G) != rank(H), the cohomology ring could have torsion.")
    
    print("\n--- Step 3: Determine the rank of the torsion subgroup ---")
    print("The rank of an abelian group is the dimension of its free part.")
    print("A torsion group consists only of elements of finite order, so its free part is trivial.")
    
    final_rank = 0
    print("\n--- Final Calculation ---")
    print(f"The rank of any torsion group is 0.")
    print(f"Therefore, the rank of the torsion subgroup is {final_rank}.")

solve_grassmannian_torsion_rank()