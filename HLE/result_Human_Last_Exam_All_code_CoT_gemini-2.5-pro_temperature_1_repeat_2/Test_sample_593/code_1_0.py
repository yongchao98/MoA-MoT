def print_treewidth_bound_formula():
    """
    This function prints the formula for the tight upper bound on the treewidth
    of the graph F, which is the union of graphs H and G intersecting on k vertices.
    """
    # Symbolic representation of the parameters
    t_H = "t_H"
    t_G = "t_G"
    k = "k"
    
    # The formula has one constant number
    number_in_equation = 1
    
    # Print the final equation for the bound
    print(f"The tight upper bound for the treewidth of F is: max({t_H}, {t_G}, {k} - {number_in_equation})")

print_treewidth_bound_formula()