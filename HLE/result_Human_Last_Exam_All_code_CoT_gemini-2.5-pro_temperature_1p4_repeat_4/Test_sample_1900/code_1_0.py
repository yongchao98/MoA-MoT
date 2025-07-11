import collections

def solve_hamiltonian_path(adj):
    """
    Solves the Hamiltonian Path problem using dynamic programming with bitmasking.
    This function both decides if a path exists and finds one if it does.

    Args:
        adj: An adjacency list representation of the graph.
             Example: [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] for K4.

    The complexity analysis is as follows:
    - n: number of vertices.
    - The dp table has 2^n * n states.
    - To compute each state, we iterate through the neighbors of a vertex.
    - The graph is sparse (max degree 8 from the problem description).
    - Let |E| be the number of edges. Total complexity is O(|E| * 2^n).
    - Since |E| = O(n), the complexity is O(n * 2^n).
    - Backtracking takes O(n) time.
    - Thus, both decision and search have a complexity of O(n * 2^n).
    """
    n = len(adj)
    if n == 0:
        print("Graph is empty.")
        return
    if n == 1:
        print("Hamiltonian path: [0]")
        return
        
    # dp[mask][u] is true if there's a path visiting vertices in 'mask' and ending at 'u'.
    dp = [[False] * n for _ in range(1 << n)]

    # Base cases: paths of length 1
    for i in range(n):
        dp[1 << i][i] = True

    # Fill the DP table for paths of length 2 to n
    for mask in range(1, 1 << n):
        for u in range(n):
            if (mask >> u) & 1:  # If u is in the current subset of vertices
                # Find a predecessor v for u in the path
                prev_mask = mask ^ (1 << u)
                if prev_mask == 0:
                    continue
                for v in adj[u]:
                    if (prev_mask >> v) & 1 and dp[prev_mask][v]:
                        dp[mask][u] = True
                        break
    
    # Check if a Hamiltonian path exists
    last_node = -1
    final_mask = (1 << n) - 1
    for i in range(n):
        if dp[final_mask][i]:
            last_node = i
            break
            
    # Reconstruct the path if one was found
    if last_node != -1:
        path = []
        current_mask = final_mask
        current_node = last_node
        
        while current_node != -1:
            path.append(current_node)
            prev_mask = current_mask ^ (1 << current_node)
            if prev_mask == 0:
                break
                
            # Find the previous node in the path
            found_prev = False
            for prev_node in adj[current_node]:
                if (prev_mask >> prev_node) & 1 and dp[prev_mask][prev_node]:
                    current_mask = prev_mask
                    current_node = prev_node
                    found_prev = True
                    break
            if not found_prev:
                 # This part should ideally not be reached if logic is correct
                 current_node = -1
        
        print("A Hamiltonian path exists.")
        print(f"Path: {path[::-1]}")
    else:
        print("No Hamiltonian path exists.")


if __name__ == '__main__':
    # Example: A 2x2 grid of red balls. This forms a K4 graph (a cycle of 4
    # with two diagonal connections). It is connected and locally connected.
    # Vertices: 0, 1, 2, 3
    # 0 -- 1
    # | \/ |
    # 2 -- 3
    # This graph has n=4 red balls.
    n_balls = 4
    # Adjacency list for K4
    adj_k4 = [
        [1, 2, 3],
        [0, 2, 3],
        [0, 1, 3],
        [0, 1, 2]
    ]
    solve_hamiltonian_path(adj_k4)
    # The final answer for the complexity questions is not computed by the code,
    # but derived from its analysis. Here we print the final string answer.
    print("\n---")
    print("Complexity Analysis Answers:")
    final_answer = "O(n * 2^n); O(n * 2^n)"
    print("Question A; Question B")
    print(final_answer)
