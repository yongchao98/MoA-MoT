import numpy as np

def calculate_total_variation():
    """
    This function demonstrates the conclusion from the problem's premises.
    We create a graph with a cycle and define vertex/edge signals that satisfy
    the given conditions. We then show that the total variation must be zero.
    """
    # Let's model a square graph with 4 vertices and 4 edges forming a cycle.
    # Vertices V = {0, 1, 2, 3}
    # Edges E = {(0,1), (1,2), (2,3), (3,0)}
    
    # From the reasoning, we concluded that the edge signal x1 must be zero.
    # x1_e = 0 for all edges e.
    x1 = np.array([0, 0, 0, 0])
    
    # The condition x1_e = |x0_u - x0_v| implies that x0 must be constant
    # across all connected vertices. Let's choose a constant value, e.g., 5.
    x0 = np.array([5, 5, 5, 5])
    
    # The edges connect vertices (0,1), (1,2), (2,3), (3,0)
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    
    # Let's verify the conditions:
    # 1. "no cycles with non-zero sum": The only cycle is the square itself.
    #    The sum of x1 values on the cycle is 0+0+0+0 = 0. Condition met.
    # 2. B1 * x1 = 0: Since x1 is the zero vector, this is trivially met.
    # 3. x1_e = |x0_u - x0_v|:
    #    |x0[0]-x0[1]| = |5-5| = 0
    #    |x0[1]-x0[2]| = |5-5| = 0
    #    |x0[2]-x0[3]| = |5-5| = 0
    #    |x0[3]-x0[0]| = |5-5| = 0
    #    This matches our x1 vector. Condition met.

    # Now, let's calculate the total variation of the graph signal x0.
    # TV = sum(|x0_u - x0_v|) over all edges {u,v}
    
    tv_terms = [abs(x0[u] - x0[v]) for u, v in edges]
    total_variation = sum(tv_terms)
    
    # Print the calculation step-by-step
    equation_str = "Total Variation = "
    term_strs = [f"|{x0[u]} - {x0[v]}|" for u, v in edges]
    equation_str += " + ".join(term_strs)
    
    value_str = " " * (len("Total Variation = "))
    value_term_strs = [f"{term}" for term in tv_terms]
    value_str += " + ".join(value_term_strs)
    
    print(equation_str)
    print(value_str)
    print(f"Total Variation = {total_variation}")

calculate_total_variation()