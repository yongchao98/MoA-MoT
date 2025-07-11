import math

def solve():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.
    """
    n = 14
    corners = {(0, 0), (0, n - 1), (n - 1, 0), (n - 1, n - 1)}
    num_corners = len(corners)

    best_config = {
        "pattern": "",
        "total_sites": float('inf'),
        "size_S": 0,
        "corners_in_S": 0
    }

    # Candidate patterns for hitting sets
    # type 0: (i+j) % 3 == k
    # type 1: (i-j) % 3 == k
    for pattern_type in range(2):
        for k in range(3):
            S = set()
            for r in range(n):
                for c in range(n):
                    if pattern_type == 0:
                        val = (r + c) % 3
                    else: # pattern_type == 1
                        val = (r - c) % 3
                    
                    if val == k:
                        S.add((r, c))

            size_S = len(S)
            corners_in_S = len(S.intersection(corners))
            
            # |I_0| = |S union C| = |S| + |C| - |S intersect C|
            total_sites = size_S + num_corners - corners_in_S
            
            if total_sites < best_config["total_sites"]:
                if pattern_type == 0:
                    p_str = f"(i+j) % 3 == {k}"
                else:
                    p_str = f"(i-j) % 3 == {k}"
                best_config["pattern"] = p_str
                best_config["total_sites"] = total_sites
                best_config["size_S"] = size_S
                best_config["corners_in_S"] = corners_in_S

    print(f"For n = {n}, the grid is {n}x{n}.")
    print("The required initial set must contain the 4 corners and hit every 2x2 square.")
    print("We search for an optimal hitting set S to minimize |S| + |Corners| - |S intersect Corners|.")
    print("\nThe best configuration found is based on the pattern:")
    print(f"S = {{(i,j) | {best_config['pattern']}}}")
    print("\nCalculating the total number of initial sites:")
    print(f"Size of the hitting set |S|: {best_config['size_S']}")
    print(f"Number of corners |C|: {num_corners}")
    print(f"Number of corners in the hitting set |S intersect C|: {best_config['corners_in_S']}")
    print(f"Total sites = |S| + |C| - |S intersect C| = {best_config['size_S']} + {num_corners} - {best_config['corners_in_S']} = {best_config['total_sites']}")
    print(f"\nThe minimum number of initially infected sites is {best_config['total_sites']}.")

solve()