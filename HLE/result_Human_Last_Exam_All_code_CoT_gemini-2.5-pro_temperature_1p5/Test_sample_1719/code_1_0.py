def get_obstruction_groups(n, k):
    """
    Prints the list of groups that classify the obstructions.

    Args:
        n: The dimension of the suspension Sigma X. X is a homology (n-1)-sphere.
        k: The bundle E is a rank 2k vector bundle.
    """
    
    # Check for valid inputs, n >= 1 for X to be a space.
    if n < 1:
        print("n must be >= 1.")
        return
    if k < 1:
        print("k must be >= 1 for SO(2k).")
        return
        
    # The groups that form the obstruction theory are determined by the
    # topology of the base space X and the fiber SO(2k).

    # The list of fundamental groups as identified in the derivation
    group_list = [
        f"H_{n-1}(X)",
        f"pi_1(SO({2*k}))",
        f"pi_{n}(SO({2*k}))"
    ]
    
    print("The homotopy-theoretic obstructions are classified by a group, which is constructed from the following list of groups:")
    for group in group_list:
        print(f"- {group}")

# Example from the problem description, where n and k are symbolic.
# To make this runnable, we can use example values. Let's say n=4, k=3.
n_val = 4
k_val = 3
get_obstruction_groups(n_val, k_val)