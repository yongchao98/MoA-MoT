import networkx as nx
from itertools import combinations

def solve_count_ans_pie(G, k):
    """
    Solves the CountAns problem using the Principle of Inclusion-Exclusion.
    This algorithm is exact but not fixed-parameter tractable.
    G: a networkx graph.
    k: a positive integer.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    total_answers = 0
    
    equation_parts = []

    # Iterate over all non-empty subsets Y of V
    for i in range(1, n + 1):
        for Y in combinations(nodes, i):
            Y = list(Y)
            
            # Compute N(Y) = intersection of neighborhoods of y in Y
            if not G.neighbors(Y[0]):
                common_neighbors = set()
            else:
                common_neighbors = set(G.neighbors(Y[0]))
            
            for j in range(1, len(Y)):
                common_neighbors.intersection_update(G.neighbors(Y[j]))
            
            num_common_neighbors = len(common_neighbors)
            term = (num_common_neighbors ** k)
            
            if term == 0:
                continue

            if len(Y) % 2 == 1:  # |Y| is odd
                total_answers += term
                equation_parts.append(f"+ {term}")
            else:  # |Y| is even
                total_answers -= term
                equation_parts.append(f"- {term}")

    print(f"Calculating for k={k} using PIE on a graph with {n} vertices.")
    # To keep the output readable, we only show the final equation summary
    final_equation = " ".join(equation_parts).lstrip("+ ")
    if not final_equation:
        final_equation = "0"
    print(f"Calculation: {final_equation} = {total_answers}")
    return total_answers

def main():
    # Example: A star graph with a center 'c' and 3 leaves 'l1', 'l2', 'l3'
    G = nx.Graph()
    G.add_edges_from([('c', 'l1'), ('c', 'l2'), ('c', 'l3')])
    k = 2

    # Manual calculation for this graph:
    # The set of answers is N('c')^k U N('l1')^k U N('l2')^k U N('l3')^k
    # = {'l1','l2','l3'}^k U {'c'}^k U {'c'}^k U {'c'}^k
    # These sets of tuples are disjoint.
    # Total count = |{'l1','l2','l3'}|^k + |{'c'}|^k = 3^k + 1^k.
    # For k=2, this is 3^2 + 1 = 9 + 1 = 10.
    
    print(f"Running on a star graph with 4 vertices.")
    num_ans = solve_count_ans_pie(G, k)
    print(f"\nThe number of answers for k={k} is: {num_ans}")
    print(f"Manual check: 3^{k} + 1 = {3**k + 1}")

if __name__ == '__main__':
    main()
