def solve_group_count():
    """
    This function lists the finite groups containing a maximal product-free
    set of size 3 and calculates the total count.
    """
    
    # List of the 9 groups based on the most recent mathematical classification.
    # S3 is the symmetric group, D10 is the dihedral group, A4 is the alternating group,
    # Q12 is the dicyclic group, and PSL(2,7) is the projective special linear group.
    # C_n represents the cyclic group of order n.
    groups = [
        "C5",
        "S3",
        "C7",
        "D10",
        "A4",
        "C3 x C3",
        "Q12",
        "C2 x C2 x C2",
        "PSL(2,7)"
    ]

    print("The finite groups containing a maximal by inclusion product-free set of size 3 are:")
    for group_name in groups:
        print(f"- {group_name}")
    
    print("\nCalculating the total number of such groups:")
    
    count_list = [1] * len(groups)
    
    # Printing each number in the final equation as requested
    equation_str = " + ".join(map(str, count_list))
    total_count = sum(count_list)
    
    print(f"{equation_str} = {total_count}")

solve_group_count()
