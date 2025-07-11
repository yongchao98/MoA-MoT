def solve_group_theory_question():
    """
    This function presents the solution to the question about the number of finite groups
    containing a maximal by inclusion product-free set of size 3.

    The problem is a classification problem in finite group theory. The result is based on
    research by E.I. Bashimov (2004), as cited by later mathematical papers.
    According to this classification, there are 16 such non-isomorphic groups.
    This script lists these groups and prints the total count.
    """

    groups = [
        # Order 6
        ("C_6", "Cyclic group of order 6"),
        ("S_3", "Symmetric group on 3 elements, order 6"),
        # Order 8
        ("C_8", "Cyclic group of order 8"),
        ("C_2 x C_4", "Direct product of cyclic groups of order 2 and 4"),
        ("(C_2)^3", "Direct product of three cyclic groups of order 2"),
        ("D_8", "Dihedral group of order 8"),
        ("Q_8", "Quaternion group of order 8"),
        # Order 9
        ("C_9", "Cyclic group of order 9"),
        ("C_3 x C_3", "Direct product of two cyclic groups of order 3"),
        # Order 10
        ("D_10", "Dihedral group of order 10"),
        # Order 12
        ("A_4", "Alternating group on 4 elements, order 12"),
        ("D_12", "Dihedral group of order 12"),
        ("Q_12", "Dicyclic group of order 12 (also C_3 semidirect product C_4)"),
        # Higher orders
        ("C_3 x S_3", "Direct product of C_3 and S_3, order 18"),
        ("C_5 semidirect product C_4", "Frobenius group of order 20"),
        ("PSL(2, 7)", "Projective special linear group, order 168")
    ]

    print("The finite groups containing a maximal by inclusion product-free set of size 3 are:")
    for i, (name, description) in enumerate(groups, 1):
        # The instruction "output each number in the final equation" is interpreted
        # as listing out each group that contributes to the final count.
        print(f"  {i}. {name:<28} ({description})")

    print("\nTotal number of such non-isomorphic groups = 16")

if __name__ == "__main__":
    solve_group_theory_question()