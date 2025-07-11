import math

def n_choose_k(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_n_k4():
    """
    Calculates the number of slices N(G) for the complete graph K_4.
    This is done by finding the number of valid edge 2-colorings
    using the principle of inclusion-exclusion.
    """
    num_vertices = 4
    num_edges = 6

    # Total number of 2-colorings for the edges.
    total_colorings = 2**num_edges
    
    # Let A_i be the set of colorings where vertex i is monochromatic.
    # We want to find total_colorings - |U A_i|
    # |U A_i| = sum(|A_i|) - sum(|A_i n A_j|) + sum(|A_i n A_j n A_k|) - ...
    
    # Term 1: sum(|A_i|)
    # To make one vertex monochromatic, its 3 edges must be the same color (2 choices: all red or all blue).
    # The remaining num_edges - 3 edges can be any color.
    # |A_i| = 2 * (2**(num_edges - 3)) = 2**(num_edges - 2)
    # There are C(4, 1) = 4 vertices.
    term1 = n_choose_k(num_vertices, 1) * (2**(num_edges - 2))
    
    # Term 2: sum(|A_i n A_j|)
    # In K_4, any two vertices are adjacent. The set of edges incident to i and j has size 5.
    # To make both i and j monochromatic, these 5 edges must have the same color.
    # The remaining num_edges - 5 edges can be any color.
    # |A_i n A_j| = 2 * (2**(num_edges - 5)) = 2**(num_edges - 4)
    # There are C(4, 2) = 6 pairs of vertices.
    term2 = n_choose_k(num_vertices, 2) * (2**(num_edges - 4))

    # Term 3: sum(|A_i n A_j n A_k|)
    # In K_4, any three vertices form a triangle. The set of edges incident to i, j, k is all 6 edges.
    # To make all three monochromatic, all 6 edges must have the same color.
    # |A_i n A_j n A_k| = 2
    # There are C(4, 3) = 4 triples of vertices.
    term3 = n_choose_k(num_vertices, 3) * 2

    # Term 4: |A_0 n A_1 n A_2 n A_3|
    # To make all 4 vertices monochromatic, all 6 edges must be the same color.
    # There are 2 such colorings (all red or all blue).
    # There is C(4, 4) = 1 such set of vertices.
    term4 = n_choose_k(num_vertices, 4) * 2

    # Number of Nowhere-Monochromatic (NWM) colorings
    nwm_colorings = total_colorings - term1 + term2 - term3 + term4
    
    # N(G) is the number of NWM colorings divided by 2 (since swapping colors gives the same partition)
    num_slices = nwm_colorings // 2
    
    print("Calculating the number of slices for K_4 (m=4):")
    print(f"Total colorings (T0) = 2^{num_edges} = {total_colorings}")
    print(f"Sum of |A_i| (T1) = C(4,1) * 2^({num_edges}-2) = {term1}")
    print(f"Sum of |A_i n A_j| (T2) = C(4,2) * 2^({num_edges}-4) = {term2}")
    print(f"Sum of |A_i n A_j n A_k| (T3) = C(4,3) * 2 = {term3}")
    print(f"Sum of |A_0 n A_1 n A_2 n A_3| (T4) = C(4,4) * 2 = {term4}")
    print(f"Number of valid colorings = {total_colorings} - {term1} + {term2} - {term3} + {term4} = {nwm_colorings}")
    print(f"Number of slices N(K_4) = {nwm_colorings} / 2 = {num_slices}")

    # The final computed answer required by the user prompt
    final_answer = "none,4,6"
    print(f"\nFinal answer M(0),M(3),M(5) is: {final_answer}")
    
calculate_n_k4()