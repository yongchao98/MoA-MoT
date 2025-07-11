def calculate_minimum_diameter(n, m):
    """
    Calculates and prints the minimum possible diameter of an undirected tree G
    with n+2 vertices and m leaves.

    Args:
        n (int): A positive integer.
        m (int): A positive integer representing the number of leaves.
    """
    print(f"--- Calculating for n={n}, m={m} ---")
    
    # Input validation based on graph theory principles.
    # 1. n and m must be positive.
    if n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        print("-" * 35 + "\n")
        return
        
    # 2. A tree with a diameter needs at least 2 leaves.
    if m < 2:
        print("Error: A tree must have at least 2 leaves to define a diameter.")
        print("-" * 35 + "\n")
        return

    # 3. For a tree with N > 2 vertices, the number of leaves is at most N-1 (for a star graph).
    # Here N = n+2. So, m <= n+1.
    if m > n + 1:
        print(f"Error: A tree with {n+2} vertices cannot have more than {n+1} leaves.")
        print("-" * 35 + "\n")
        return

    # The minimum possible diameter is given by the formula: min(4, n - m + 3)
    # This formula is derived from analyzing the structure of the tree's internal nodes.

    # Step 1: Calculate the value of the expression `n - m + 3`
    expression_value = n - m + 3
    
    # Step 2: The diameter is the minimum of 4 and the value from Step 1.
    diameter = min(4, expression_value)
    
    # Printing the steps of the calculation as requested
    print(f"The minimum diameter 'd' is calculated using the formula: d = min(4, n - m + 3)")
    print(f"Substituting the values of n and m:")
    print(f"d = min(4, {n} - {m} + 3)")
    print(f"d = min(4, {expression_value})")
    print(f"The final minimum possible diameter is: {diameter}")
    print("-" * 35 + "\n")

# --- Main execution ---
# Running the function for three different cases to demonstrate the logic.

# Case 1: m = n + 1 (e.g., n=8, m=9) -> should result in diameter 2
calculate_minimum_diameter(8, 9)

# Case 2: m = n (e.g., n=8, m=8) -> should result in diameter 3
calculate_minimum_diameter(8, 8)

# Case 3: m <= n - 1 (e.g., n=8, m=5) -> should result in diameter 4
calculate_minimum_diameter(8, 5)

# Another example for Case 3 to show it's always 4
calculate_minimum_diameter(10, 3)
