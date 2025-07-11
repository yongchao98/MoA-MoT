import math

def compute_l2_betti_number():
    """
    This function computes the first l^2-Betti number of the group G
    based on the Chiswell-Lück formula for graphs of groups.
    """

    # --- 1. Graph Term ---
    # The underlying graph Y is the line graph of the Petersen graph.
    # Petersen graph has 10 vertices (all degree 3) and 15 edges.
    # Line graph Y = L(P) has:
    V_Y = 15
    E_Y = 10 * math.comb(3, 2) # Each vertex in P has degree 3
    
    # The fundamental group pi_1(Y) is a free group F_r.
    # The rank r is the first Betti number of the graph Y.
    rank_pi1_Y = E_Y - V_Y + 1
    
    # For a free group F_r (r>=1), beta_0 = 0 and beta_1 = r - 1.
    beta0_pi1_Y = 0
    beta1_pi1_Y = rank_pi1_Y - 1
    
    graph_term = beta1_pi1_Y - beta0_pi1_Y

    # --- 2. Edge Group Term ---
    # Each edge group G_e is isomorphic to some N_g.
    # N_g is a finite-index subgroup of the fundamental group of a hyperbolic 3-manifold.
    # For such groups, all l^2-Betti numbers are 0.
    # beta_k(N_g) = 0 for all k, g.
    # N_g is infinite, so beta_0(N_g) = 0.
    # Thus, the contribution for each edge is beta_1(G_e) - beta_0(G_e) = 0 - 0 = 0.
    sum_e_term = 0

    # --- 3. Vertex Group Term ---
    # This is the sum of (beta_1(G_v) - beta_0(G_v)) over all 15 vertices.
    sum_v_term = 0
    
    # For v_1, G_{v_1} = N_{100}. The term is beta_1(N_100) - beta_0(N_100) = 0 - 0 = 0.
    term_v1 = 0
    sum_v_term += term_v1

    # For v_i where i = 2, ..., 15:
    # G_{v_i} is the free product of k = i-1 groups (N_2, ..., N_i).
    # The formula for a free product is:
    # beta_1(*) - beta_0(*) = sum(beta_1(H_j) - beta_0(H_j)) + k - 1
    # For G_{v_i}, beta_0(G_{v_i}) = 0 as it's an infinite group.
    # And beta_1(N_g) = 0, beta_0(N_g) = 0.
    # So beta_1(G_{v_i}) = sum(0-0) + (i-1) - 1 = i - 2.
    # The contribution is beta_1(G_{v_i}) - beta_0(G_{v_i}) = (i-2) - 0 = i-2.
    sum_v_others = sum(i - 2 for i in range(2, 16))
    sum_v_term += sum_v_others

    # --- 4. Final Calculation ---
    # The group G is infinite, so beta_0(G) = 0. The formula is:
    # beta_1(G) = sum_v_term - sum_e_term + graph_term
    result = sum_v_term - sum_e_term + graph_term

    print("Step-by-step calculation for the first l^2-Betti number of G:")
    print("-" * 60)
    print("The Chiswell-Lück formula is:")
    print("β₁⁽²⁾(G) - β₀⁽²⁾(G) = Σᵥ(β₁⁽²⁾(Gᵥ) - β₀⁽²⁾(Gᵥ)) - Σₑ(β₁⁽²⁾(Gₑ) - β₀⁽²⁾(Gₑ)) + (β₁⁽²⁾(π₁Y) - β₀⁽²⁾(π₁Y))")
    print("\nSince G is infinite, β₀⁽²⁾(G) = 0. The formula simplifies to:")
    print("β₁⁽²⁾(G) = [Vertex Sum] - [Edge Sum] + [Graph Term]")
    print("-" * 60)
    
    print(f"1. The Vertex Sum is calculated as Σᵢ(i-2) for i from 2 to 15:")
    print(f"   Σᵥ(β₁⁽²⁾(Gᵥ) - β₀⁽²⁾(Gᵥ)) = {sum_v_term}")

    print(f"\n2. The Edge Sum is zero as all βₖ⁽²⁾(N_g) are zero:")
    print(f"   Σₑ(β₁⁽²⁾(Gₑ) - β₀⁽²⁾(Gₑ)) = {sum_e_term}")

    print(f"\n3. The Graph Term is calculated from the graph Y = L(Petersen):")
    print(f"   |V(Y)|={V_Y}, |E(Y)|={E_Y}, rank(π₁Y)={rank_pi1_Y}")
    print(f"   β₁⁽²⁾(π₁Y) - β₀⁽²⁾(π₁Y) = {beta1_pi1_Y} - {beta0_pi1_Y} = {graph_term}")

    print("-" * 60)
    print("Putting it all together:")
    print(f"β₁⁽²⁾(G) = {sum_v_term} - {sum_e_term} + {graph_term} = {result}")

compute_l2_betti_number()