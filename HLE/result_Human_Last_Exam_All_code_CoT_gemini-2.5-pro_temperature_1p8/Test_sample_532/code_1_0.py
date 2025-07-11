def classify_finite_filled_nilpotent_groups():
    """
    This function explains and presents the classification of finite filled nilpotent groups
    based on established mathematical research.
    """

    print("--- The Classification of Finite Filled Nilpotent Groups ---\n")
    print("A finite group G is nilpotent if it is the direct product of its Sylow p-subgroups.")
    print("A finite nilpotent group is 'filled' if and only if each of its Sylow p-subgroups is filled.\n")
    print("Based on research in group theory, a p-group can only be filled if p is 2, 3, or 5.")
    print("Therefore, a finite nilpotent group G is filled if and only if it is a direct product of the form:")
    print("G = G_2 x G_3 x G_5, where G_p is a filled Sylow p-subgroup (which can be the trivial group).\n")
    print("The possible structures for these filled p-subgroups are as follows:\n")

    filled_2_groups = [
        "C_2 (Cyclic group of order 2)",
        "C_4 (Cyclic group of order 4)",
        "C_2 x C_2 (Klein four-group)",
        "C_4 x C_2 (Abelian group of order 8)",
        "D_8 (Dihedral group of order 8)",
        "Q_8 (Quaternion group of order 8)"
    ]

    filled_3_groups = [
        "C_3 (Cyclic group of order 3)",
        "C_3 x C_3 (Abelian group of order 9 and exponent 3)"
    ]

    filled_5_groups = [
        "C_5 (Cyclic group of order 5)",
        "C_5 x C_5 (Abelian group of order 25 and exponent 5)"
    ]

    print("1. For primes p >= 7:")
    print("   - The Sylow p-subgroup must be the trivial group {e}.\n")

    print("2. For p = 2, the Sylow 2-subgroup (G_2) must be trivial or one of:")
    for group in filled_2_groups:
        print(f"   - {group}")
    print()

    print("3. For p = 3, the Sylow 3-subgroup (G_3) must be trivial or one of:")
    for group in filled_3_groups:
        print(f"   - {group}")
    print()

    print("4. For p = 5, the Sylow 5-subgroup (G_5) must be trivial or one of:")
    for group in filled_5_groups:
        print(f"   - {group}")
    print()

    print("Any direct product of one choice from each of these categories (including the trivial choice) results in a finite filled nilpotent group.")


if __name__ == '__main__':
    classify_finite_filled_nilpotent_groups()
