def display_treewidth_bound_formula():
    """
    Displays the formula for the tight upper bound on the treewidth of the union of two graphs.
    """
    # The variables are symbolic, as no numerical values were provided in the problem.
    t_H = "t_H"
    t_G = "t_G"
    k = "k"
    one = "1"

    print("The tight upper bound on the treewidth of the graph F, denoted t_F, is given by:")
    
    # We output each term of the equation as requested.
    print(f"t_F <= max({t_H}, {t_G}, {k} - {one})")

if __name__ == '__main__':
    display_treewidth_bound_formula()