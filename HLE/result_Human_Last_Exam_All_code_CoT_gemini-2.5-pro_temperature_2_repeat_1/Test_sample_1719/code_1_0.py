def print_obstruction_groups():
    """
    This function prints the list of homology and homotopy groups that form the obstructions
    for the two paths of bundle automorphisms to be homotopic.
    """
    
    n_str = "n"
    k_str = "k"
    
    print("The obstructions are classified by the group pi_1(Aut(E)), which is constructed from the following groups:")
    
    groups = [
        f"The reduced homology group H_{{{n_str}-1}}(X), which is Z.",
        f"pi_1(SO(2*{k_str}))",
        f"pi_2(SO(2*{k_str}))",
        f"pi_{n_str}(SO(2*{k_str}))",
        f"pi_{n_str}+1(SO(2*{k_str}))"
    ]
    
    for i, group in enumerate(groups):
        # We replace the pythonic f-string formatting with mathematical notation for clarity.
        group_str = group.replace("pi_", "π_").replace("*", "")
        group_str = group_str.replace("{", "{").replace("}", "}")
        group_str = group_str.replace("H_{n-1}(X), which is Z", "H̃_{n-1}(X)")
        print(f"  {i+1}. {group_str}")

print_obstruction_groups()