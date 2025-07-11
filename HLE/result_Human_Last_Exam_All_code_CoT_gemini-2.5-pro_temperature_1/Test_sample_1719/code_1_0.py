def print_obstruction_groups(n, k):
    """
    Prints the list of homotopy-theoretic obstruction groups.

    Args:
        n: The dimension of the sphere S^(n-1) that X is a homology sphere of.
        k: Half the rank of the vector bundle E (rank is 2k).
    """
    print("The homotopy-theoretic obstructions are elements of the following list of groups:")
    print("-" * 70)

    # First obstruction from the fibration
    print(f"1. pi_1(SO({2*k}))")

    # Higher obstructions from obstruction theory
    print("\n2. The sequence of higher-order obstructions, which lie in the groups:")
    # p runs from 3 to n+1
    for p in range(3, n + 1 + 1):
        # The obstruction group is H^p(Sigma^2 X, pi_p(SO(2k)))
        # which is isomorphic to H_tilde^{p-2}(X, pi_p(SO(2k)))
        print(f"   ~H^{p-2}(X; pi_{p}(SO({2*k})))")

    print("-" * 70)
    print("\nNote: For a homology (n-1)-sphere X, most of these higher-order groups vanish.")
    print("Specifically, since H_tilde_i(X) is non-zero only for i=n-1, the only potentially")
    print(f"non-zero group in the second list is the last one: ~H^{n-1}(X; pi_{n+1}(SO({2*k}))).")


if __name__ == '__main__':
    # Example values for a homology 3-sphere (n-1=3 => n=4) and a rank 8 bundle (2k=8 => k=4)
    n_val = 4
    k_val = 4
    print_obstruction_groups(n_val, k_val)
    # The final answer is the list of groups printed by the function.
    # The <<<...>>> format is not applicable here as the answer is a descriptive list.
    # To conform to the requested output format, we will output the list itself as a string.
    
    obstruction_list = []
    obstruction_list.append(f"pi_1(SO({2*k_val}))")
    for p in range(3, n_val + 1 + 1):
        obstruction_list.append(f"~H^{p-2}(X; pi_{p}(SO({2*k_val})))")
    
    final_answer = str(obstruction_list)
    # print(f"\n<<< {final_answer} >>>")
    # The format seems to expect a single value. Since the problem asks for a list
    # and has a complex nature, presenting the code that generates the list is the
    # most accurate and helpful response. The final list for n=4, k=4 is:
    # ['pi_1(SO(8))', '~H^1(X; pi_3(SO(8)))', '~H^2(X; pi_4(SO(8)))', '~H^3(X; pi_5(SO(8)))']
    
    # As the final output requested is a single value, and my analysis shows the obstructions
    # are composed of two main parts, I will provide the list of the two key components that
    # are always present in the obstruction analysis.
    final_obstruction_components = [f"pi_1(SO({2*k_val}))", f"H^{n_val+1}(Sigma^2(X); pi_{n_val+1}(SO({2*k_val})))"]
    # For n=4, k=4, this would be ['pi_1(SO(8))', 'H^5(Sigma^2(X); pi_5(SO(8)))']
    
    # Let's provide the full list as the final answer in the requested format.
    print(f"\n<<< {str(obstruction_list)} >>>")