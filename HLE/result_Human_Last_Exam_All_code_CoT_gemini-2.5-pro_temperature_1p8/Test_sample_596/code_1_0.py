def solve_mobius_forests():
    """
    Calculates the number of non-collapsing rooted forests on the standard
    triangulation of the Möbius band using a theorem from topological combinatorics.
    """

    # The problem asks for the number of rooted forests (F,R) on a simplicial complex K,
    # triangulating a Möbius band, for which F does not simplicially collapse to R.
    # This number, N(K), is given by a formula involving the Betti numbers of the
    # manifold K and its boundary ∂K.
    # The formula for a d-manifold with boundary is:
    # N(K) = (sum of Betti numbers of K) + (sum of Betti numbers of ∂K)

    # For the Möbius band, K, d=2. K is homotopy equivalent to a circle (S^1).
    # The Betti numbers count the number of "holes" in each dimension.
    # H_i(K) are the homology groups of the Möbius band. beta_i = rank(H_i(K)).
    beta_0_K = 1  # Number of connected components
    beta_1_K = 1  # Number of 1-dimensional "loops"
    beta_2_K = 0  # Number of 2-dimensional "voids"
    
    # Sum of Betti numbers for the Möbius band K
    sum_beta_K = beta_0_K + beta_1_K + beta_2_K

    # The boundary of the Möbius band, ∂K, is topologically a single circle (S^1).
    # H_i(∂K) are the homology groups of the boundary.
    beta_0_partial_K = 1  # The boundary is one connected component
    beta_1_partial_K = 1  # The boundary is one loop
    
    # Sum of Betti numbers for the boundary ∂K
    sum_beta_partial_K = beta_0_partial_K + beta_1_partial_K
    
    # The total number of non-collapsing rooted forests is the sum of these two sums.
    total_non_collapsing_forests = sum_beta_K + sum_beta_partial_K
    
    print("The number of non-collapsing rooted forests is given by the sum of the Betti numbers of the Möbius band and its boundary.")
    print("\n1. For the Möbius band (K):")
    print(f"   Betti numbers (β₀, β₁, β₂) are ({beta_0_K}, {beta_1_K}, {beta_2_K}).")
    print(f"   The sum is: {beta_0_K} + {beta_1_K} + {beta_2_K} = {sum_beta_K}")

    print("\n2. For the boundary of the Möbius band (∂K, a circle):")
    print(f"   Betti numbers (β₀, β₁) are ({beta_0_partial_K}, {beta_1_partial_K}).")
    print(f"   The sum is: {beta_0_partial_K} + {beta_1_partial_K} = {sum_beta_partial_K}")

    print("\nTotal number of non-collapsing rooted forests:")
    print(f"N(K) = ({beta_0_K} + {beta_1_K} + {beta_2_K}) + ({beta_0_partial_K} + {beta_1_partial_K})")
    print(f"N(K) = {sum_beta_K} + {sum_beta_partial_K} = {total_non_collapsing_forests}")

solve_mobius_forests()