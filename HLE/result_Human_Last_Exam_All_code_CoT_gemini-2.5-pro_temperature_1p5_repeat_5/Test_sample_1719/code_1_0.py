def print_obstruction_groups(n, k):
    """
    Prints the names of the groups that classify the homotopy-theoretic
    obstructions for the given problem.

    Args:
        n (int): The dimension parameter from the (n-1)-homology sphere X.
        k (int): The parameter from the rank 2k vector bundle.
    """
    print("The homotopy-theoretic obstructions for the two paths to be homotopic are classified by a group,")
    print("which is an extension of subgroups/quotients of the following groups:")
    print(f"1. The fundamental group of the special orthogonal group: pi_1(SO({2*k}))")
    print(f"2. A higher homotopy group of the special orthogonal group, determined by the dimension of X: pi_{n+1}(SO({2*k})) which is pi_{{{n}+1}}(SO({2*k}))")
    print("\nNote: The precise structure of the obstruction group also depends on pi_2 and pi_{n-1} of SO(2k) and the specific bundle E.")

if __name__ == '__main__':
    # These are parameters from the problem statement.
    # We use placeholder values as n and k were not specified numerically.
    n_placeholder = 'n'
    k_placeholder = 'k'
    
    # We format the output strings with these placeholders.
    print("Let n be the dimension parameter from the homology (n-1)-sphere X.")
    print("Let 2k be the rank of the vector bundle E.")
    
    obstruction_groups = [
        f"pi_1(SO(2k))",
        f"pi_{n+1}(SO(2k))"
    ]
    
    print("\nThe obstructions lie in the following list of groups:")
    for group in obstruction_groups:
        print(f"- {group}")
        
    print("\nMore specifically, the complete list of groups that determine the extension problem for the obstruction group is:")
    
    full_list = [
        f"H_{n-1}(X)",
        f"pi_1(SO(2k))",
        f"pi_2(SO(2k))",
        f"pi_{n-1}(SO(2k))",
        f"pi_{n+1}(SO(2k))"
    ]
    for group in full_list:
        print(f"- {group}")
