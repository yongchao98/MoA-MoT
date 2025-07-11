def solve():
    """
    This function provides the formula for the tight upper bound on the treewidth of F.
    Since the values of t_H, t_G, and k are symbolic, we print the formula.
    """
    
    t_H = "t_H"
    t_G = "t_G"
    k = "k"
    
    # The derivation shows the tight upper bound is max(t_H + k - 1, t_G + k - 1)
    
    term1_val = f"{t_H} + {k} - 1"
    term2_val = f"{t_G} + {k} - 1"
    
    bound_expression = f"max({term1_val}, {term2_val})"
    
    print("The tight upper bound on the treewidth of F is given by the expression:")
    print(bound_expression)

solve()