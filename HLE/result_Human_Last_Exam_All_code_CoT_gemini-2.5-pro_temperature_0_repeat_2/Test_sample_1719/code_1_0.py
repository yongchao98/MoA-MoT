def get_obstruction_groups(n_str="n", k_str="k"):
    """
    This function generates the list of groups where the obstructions lie.
    The obstructions determine if two paths of bundle automorphisms are homotopic.

    Args:
        n_str (str): A string representing the dimension 'n' from the homology of X.
        k_str (str): A string representing the parameter 'k' from the bundle rank 2k.
        
    Returns:
        list: A list of strings, where each string describes an obstruction group.
    """
    
    # The primary obstruction is local and lies in the fundamental group of SO(2k).
    # The numbers/variables involved are 1 and 2*k.
    group1 = f"pi_{1}(SO(2*{k_str}))"
    
    # The secondary obstruction is global and relates the topology of the base space
    # to the homotopy of the fiber.
    # The numbers/variables involved are n-1 and n+1.
    group2 = f"Hom(H_{n_str}-1(X), pi_{n_str}+1(SO(2*{k_str})))"
    
    return [group1, group2]

# --- Main execution ---
if __name__ == "__main__":
    # Get the list of obstruction groups using symbolic variables n and k.
    obstruction_groups = get_obstruction_groups()
    
    print("The homotopy-theoretic obstructions for the two paths to be homotopic are given by elements in the following groups:")
    
    # Print each group from the list.
    # The first must vanish to be able to measure the second.
    for i, group in enumerate(obstruction_groups):
        print(f"{i+1}. {group}")
        
    # The prompt asks to output each number in the final equation.
    # Here is a breakdown of the indices and parameters in the equations above.
    print("\n--- Analysis of the numbers and indices in the equations ---")
    print("In the first group, pi_1(SO(2*k)):")
    print(" - The homotopy group index is 1.")
    print(" - The rank of the special orthogonal group is 2*k.")
    
    print("\nIn the second group, Hom(H_{n-1}(X), pi_{n+1}(SO(2*k))):")
    print(" - The homology group index is n-1.")
    print(" - The homotopy group index is n+1.")
    print(" - The rank of the special orthogonal group is 2*k.")
