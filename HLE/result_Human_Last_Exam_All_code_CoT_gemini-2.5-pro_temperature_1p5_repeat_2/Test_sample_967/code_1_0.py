import math

def solve_betti_number():
    """
    Calculates the first l^2-Betti number of the fundamental group G.
    """
    
    # According to the plan:
    # 1. Calculate the sum of the l^2-Betti numbers of the vertex groups.
    # The contribution from v_1 is 0.
    # The contribution from v_i for i in [2, 15] is (i-2).
    
    # We sum (i-2) for i from 2 to 15. This is equivalent to summing j from 0 to 13.
    n = 13
    sum_vertex_betti = n * (n + 1) // 2
    
    # 2. Calculate the sum of the l^2-Betti numbers of the edge groups.
    # As explained, the first l^2-Betti number of each edge group N_g is 0.
    # Since the sum of zeros is zero, this total is 0.
    sum_edge_betti = 0
    
    # 3. Apply the Chiswell-Lück formula.
    final_betti_number = sum_vertex_betti - sum_edge_betti
    
    # Print the components of the final equation as requested.
    print(f"Sum of vertex group l^2-Betti numbers (Σ β₁(Gᵥ)): {sum_vertex_betti}")
    print(f"Sum of edge group l^2-Betti numbers (Σ β₁(Gₑ)): {sum_edge_betti}")
    print(f"The final equation is: {sum_vertex_betti} - {sum_edge_betti} = {final_betti_number}")
    print(f"The computed first l^2-Betti number of G is: {final_betti_number}")

solve_betti_number()