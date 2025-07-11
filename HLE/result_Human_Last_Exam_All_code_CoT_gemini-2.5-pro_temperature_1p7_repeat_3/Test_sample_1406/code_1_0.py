import collections
import itertools

def has_cut_points(n, k=4):
    """
    Checks if the n-cube has cut points using a grid approximation.
    n: dimension of the cube.
    k: number of divisions per axis (k+1 points). k>=3 is needed for interior points.
    """
    if n == 0:
        return False # A single point has no cut points.

    # Generate vertices of the grid graph for the n-cube
    # Vertices are n-tuples of integers from 0 to k.
    nodes = list(itertools.product(range(k + 1), repeat=n))
    num_total_nodes = len(nodes)

    # An interior point p=(p_1,...,p_n) is one where 0 < p_i < k for all i.
    # We only need to check one interior point due to symmetry.
    if k < 2: # No interior points if k < 2
        return False
        
    p_to_remove = tuple(k // 2 for _ in range(n))
    
    # We need to check if the graph without p_to_remove is connected.
    # We can do this with a graph traversal like Breadth-First Search (BFS).
    
    # Select a starting node for the BFS, ensuring it's not the removed one.
    start_node = tuple(0 for _ in range(n))
    if start_node == p_to_remove: # Should not happen with our choice of p
        start_node = tuple(1 for _ in range(n))

    q = collections.deque([start_node])
    visited = {start_node}
    
    while q:
        current_node = q.popleft()
        
        # Find neighbors of the current node
        for i in range(n):
            for move in [-1, 1]:
                neighbor_list = list(current_node)
                neighbor_list[i] += move
                neighbor = tuple(neighbor_list)
                
                # Check if the neighbor is valid (within the grid)
                if all(0 <= c <= k for c in neighbor):
                    if neighbor != p_to_remove and neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)
                        
    # If the number of visited nodes is less than the total number of nodes
    # minus the one we removed, the graph is disconnected.
    return len(visited) < num_total_nodes - 1

def main():
    """
    Main function to analyze n-cubes for n=1, 2, 3, 4.
    """
    print("Analyzing n-cubes for the cut-point property.")
    print("A space has cut points if removing a single point can disconnect it.")
    print("-" * 30)

    failing_n = []
    for n in range(1, 5):
        if has_cut_points(n):
            print(f"For n={n}, the n-cube [0,1]^{n} HAS cut points.")
            failing_n.append(n)
        else:
            print(f"For n={n}, the n-cube [0,1]^{n} does NOT have cut points.")
    
    print("-" * 30)
    print("A theorem states that a set of non-block points cannot have cut points.")
    print("Therefore, the n-cube fails to be a set of non-block points only for the values of n listed above.")
    
    count = len(failing_n)
    print(f"\nThe n-cube fails to occur as the set of non-block points for {count} value(s) of n.")
    
    # Final equation part
    print("\nThe problem is to find the total count of such values of n.")
    print(f"The set of failing n is: {failing_n}")
    final_equation = f"Total Count = {count}"
    print(final_equation)
    # Printing each number in the final equation:
    print(f"The number in the final result is: {count}")


if __name__ == "__main__":
    main()
