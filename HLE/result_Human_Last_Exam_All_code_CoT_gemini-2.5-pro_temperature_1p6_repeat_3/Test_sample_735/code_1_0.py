def solve_group_problem():
    """
    Solves the problem by presenting the known classification of finite groups
    with maximal product-free sets of size 3.
    """
    
    # According to published research in combinatorial group theory, the following
    # is the complete list of non-isomorphic finite groups containing a
    # maximal by inclusion product-free set of size 3.
    groups = [
        "C5, the cyclic group of order 5",
        "C7, the cyclic group of order 7",
        "S3, the symmetric group on 3 elements (isomorphic to D6, the dihedral group of order 6)",
        "C2 x C4, the direct product of cyclic groups of order 2 and 4",
        "Q8, the quaternion group of order 8",
        "D8, the dihedral group of order 8"
    ]
    
    print("The finite groups containing maximal by inclusion product-free sets of size 3 are:")
    for group_name in groups:
        print(f"- {group_name}")
    
    print("\nThe total number of such groups is:")
    print(len(groups))

if __name__ == "__main__":
    solve_group_problem()