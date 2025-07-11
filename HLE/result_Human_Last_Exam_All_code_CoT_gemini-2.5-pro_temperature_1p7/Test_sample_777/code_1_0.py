import networkx as nx
import math

def find_disjoint_cycles_in_grid(width, height, k, min_len):
    """
    Finds k vertex-disjoint cycles in a grid graph G.

    This function demonstrates that a sufficiently large grid contains k disjoint
    cycles of a certain minimum length. It's part of the argument for why
    the DisjointCycles problem is FPT.

    Args:
        width (int): The width of the grid.
        height (int): The height of the grid.
        k (int): The number of vertex-disjoint cycles to find.
        min_len (int): The minimum length of each cycle.

    Returns:
        A list of lists, where each inner list contains the vertices of a
        found cycle. Returns an empty list if the grid is too small.
    """
    print(f"Attempting to find {k} disjoint cycles of length >= {min_len} in a {width}x{height} grid.")

    # Check if the grid is large enough for our construction method.
    # We need 2k rows to make k disjoint cycles using pairs of rows.
    # We need 2*width to be the cycle length, so 2*width >= min_len.
    if height < 2 * k or width < math.ceil(min_len / 2.0):
        print("Grid is too small for this construction method.")
        return []

    G = nx.grid_2d_graph(width, height)
    found_cycles = []
    
    print("Grid is large enough. Constructing cycles...")
    for i in range(k):
        cycle = []
        # We use rows 2*i and 2*i+1 for the i-th cycle.
        row1 = 2 * i
        row2 = 2 * i + 1

        # Go right along the top row of the pair
        for x in range(width):
            cycle.append((x, row1))
        
        # Go left along the bottom row of the pair
        for x in range(width - 1, -1, -1):
            cycle.append((x, row2))
            
        found_cycles.append(cycle)

    return found_cycles

if __name__ == '__main__':
    # Problem instance: k=3. We need 3 disjoint cycles of length at least 3.
    k_param = 3
    min_len_param = 3 # This is k in the problem statement
    
    # Let's determine the required grid size for our method
    # Height should be at least 2*k = 6
    # Width should be at least ceil(k/2) = ceil(3/2) = 2
    grid_w = 4 # Let's make it a bit larger than minimal
    grid_h = 6

    cycles = find_disjoint_cycles_in_grid(grid_w, grid_h, k_param, min_len_param)

    if cycles:
        print(f"\nSuccessfully found {len(cycles)} disjoint cycles.")
        for i, cycle_nodes in enumerate(cycles):
            # The prompt asked to "output each number in the final equation"
            # We interpret this as printing the nodes of each found cycle.
            print(f"Cycle {i+1} (length {len(cycle_nodes)}):")
            print(' -> '.join(map(str, cycle_nodes)) + ' -> ' + str(cycle_nodes[0]))
            
            # Verify they are vertex-disjoint by converting to sets
            if i > 0:
                set1 = set(cycles[i-1])
                set2 = set(cycles[i])
                assert len(set1.intersection(set2)) == 0, "Cycles are not disjoint!"
    else:
        print("Could not find the cycles with the given parameters.")