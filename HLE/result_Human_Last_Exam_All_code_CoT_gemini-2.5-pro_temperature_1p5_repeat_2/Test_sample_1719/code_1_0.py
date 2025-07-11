def get_obstruction_groups(n, k):
    """
    Prints the list of groups that classify the homotopy-theoretic obstructions.

    Args:
        n (int): The dimension of the suspension, n = dim(X) + 1.
        k (int): Half the rank of the vector bundle E.
    """
    print("The obstructions for the two paths to be homotopic are classified by elements in the following list of groups:")
    # First obstruction group
    print(f"1. pi_1(SO({2*k}))")
    # Higher obstructions from the non-triviality of the map over SigmaX
    print(f"2. pi_{n-1}(SO({2*k}))")
    print(f"3. pi_{n}(SO({2*k}))")

# Example for a homology 3-sphere (n=4) and a rank 8 bundle (k=4)
# n = (dim X) + 1, so for a homology 3-sphere X, n = 3 + 1 = 4.
n_example = 4 
k_example = 4
get_obstruction_groups(n_example, k_example)
