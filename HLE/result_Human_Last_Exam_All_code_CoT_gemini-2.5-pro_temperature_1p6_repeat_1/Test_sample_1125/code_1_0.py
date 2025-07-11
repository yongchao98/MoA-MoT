def solve_hat_puzzle():
    """
    This script explains the solution to the 12-person hat puzzle
    by applying graph theory concepts.
    """
    
    # 1. Define the parameters of the puzzle
    num_members = 12
    
    print("Solving the 12-person hat puzzle...")
    print(f"There are {num_members} team members and {num_members} unique hats.")
    print("The goal is to find N, the maximum number of people guaranteed to know their hat number.\n")

    # 2. Explain the optimal strategy using a graph model
    print("--- The Strategy ---")
    print("1. The members are vertices in a graph, and pairs who raise hands form edges.")
    print("2. To guarantee information can be shared for deduction, the graph must be connected.")
    print(f"3. The team makes {num_members - 1} pairing requests to form a connected graph (a tree).\n")
    
    # 3. Explain the leader's role and the concept of a vertex cover
    print("--- The Leader's Counter-Strategy ---")
    print("1. When a pair (A,B) asks, the leader reveals the number of A or B.")
    print("2. The set of people whose numbers are revealed must 'cover' every pair.")
    print("   This is a 'vertex cover' in graph theory.")
    print("3. To minimize the number of winners, the leader will choose revelations that form a")
    print("   'minimum vertex cover' (the smallest possible cover), of size τ(G).\n")

    # 4. Explain how the team maximizes their guaranteed successes
    print("--- Finding the Optimal N ---")
    print("1. The team must choose a tree graph (G) that maximizes the size of the minimum vertex cover, τ(G).")
    print(f"2. We use the identity: τ(G) + α(G) = n, where n = {num_members} and α(G) is the 'independence number'.")
    print("   This means maximizing τ(G) is the same as minimizing α(G).")
    print(f"3. For any tree on n vertices, it is known that α(G) ≥ n / 2.")
    
    # 5. Calculate and present the final equation and answer
    min_alpha = num_members // 2
    max_tau = num_members - min_alpha
    
    print("\n--- The Final Calculation ---")
    print("The final equation is derived from these principles:")
    print(f"Minimum α(G) for a tree on {num_members} vertices = {min_alpha}")
    print(f"Maximum τ(G) = {num_members} - Minimum α(G)")
    print(f"Max N = {num_members} - {min_alpha} = {max_tau}")

    print("\nThis maximum is achievable if the team forms a simple 'path' or 'line' graph.")
    print("The leader can then guarantee revealing the numbers of 6 people, but no fewer.")
    print("The other 6 will not have their numbers revealed and won't have enough information to deduce them.")
    print("\nTherefore, the largest possible value of N is 6.")

# Run the explanation
solve_hat_puzzle()