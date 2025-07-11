def get_obstruction_groups():
    """
    This function prints the list of homotopy groups where the obstructions lie.
    The parameters n and k from the problem are treated symbolically as they are not given specific values.
    """
    # The rank of the special orthogonal group is 2k. This is symbolic.
    rank_str = "2k"

    # The first obstruction group is the fundamental group of SO(2k). Its index is 1.
    index1 = 1

    # The second obstruction group's index depends on n, the dimension parameter of the homology sphere.
    # The index is n+1. We will print the numbers involved.
    index2_n_part_str = "n"
    index2_const_part = 1
    
    print("The homotopy obstructions lie in a group that is the direct sum of the following two groups:")
    print("-" * 70)
    # Print the first group
    print(f"1. The fundamental group of SO({rank_str}): pi_{index1}(SO({rank_str}))")
    
    # Print the second group
    print(f"2. A higher homotopy group of SO({rank_str}): pi_{{{index2_n_part_str} + {str(index2_const_part)}}}(SO({rank_str}))")
    print("-" * 70)
    print("The index 'n' is such that X is a homology (n-1)-sphere.")

# Execute the function to display the result.
get_obstruction_groups()