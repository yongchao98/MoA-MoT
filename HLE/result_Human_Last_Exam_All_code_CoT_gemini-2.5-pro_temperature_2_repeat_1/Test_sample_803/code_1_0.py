def solve_filled_groups_problem():
    """
    This function outlines the solution to find all nonabelian filled groups
    of order 2*q^m for an odd prime q and natural number m.
    The solution is based on established theorems in group theory.
    """

    print("--- Analysis of Nonabelian Filled Groups of Order 2*q^m ---")
    print("\nIntroduction:")
    print("A group G is called 'filled' if every maximal product-free set in G fills G.")
    print("We want to find all nonabelian filled groups of order 2*q^m, where q is an odd prime and m is a natural number.")
    print("The solution relies on a classification theorem for non-abelian filled groups of even order.")

    print("\nKey Theorem:")
    print("A non-abelian group G of even order is filled if and only if it has no quotient group isomorphic to any of the following 'forbidden' groups:")
    forbidden_groups = [
        "C_{3k} for any integer k >= 1",
        "C_2 x C_4 (order 8)",
        "D_10, the dihedral group of order 10",
        "D_6 (isomorphic to S_3), the dihedral group of order 6",
        "Dic_12, the dicyclic group of order 12",
        "C_3 x C_3 (order 9)",
        "S_3 x C_3 (order 18)",
        "(C_3 x C_3) semi-direct product C_2 (order 18)",
        "A_4, the alternating group of order 12"
    ]
    print(" - " + "\n - ".join(forbidden_groups))
    print("\nThe order of a quotient of G must divide the order of G. So, we check which of these forbidden groups can be a quotient of a group of order 2*q^m.")
    print("The prime factors of the orders of the forbidden groups are in the set {2, 3, 5}.")

    print("\n--- Case-by-Case Analysis Based on q ---")

    # Case 1: q = 3
    print("\nCase 1: q = 3")
    print("The order of the group is 2 * 3^m.")
    print("A group of this order can potentially have quotients with prime factors 2 and 3, such as S_3 or C_{3k}.")
    print("Let G be a non-abelian group of order 2 * 3^m. G has a normal Sylow 3-subgroup, say Q.")
    print("It can be shown that G must have a quotient group isomorphic to either S_3 or an abelian group whose order is divisible by 3 (like C_3 or C_6). Both S_3 and C_{3k} are on the forbidden list.")
    print("Therefore, for q = 3, there are NO nonabelian filled groups of order 2 * 3^m.")

    # Case 2: q = 5
    print("\nCase 2: q = 5")
    print("The order of the group is 2 * 5^m.")
    print("Checking the forbidden list, the only group whose order has prime factors from {2, 5} is D_10 (order 10).")
    print("So, a group G of order 2 * 5^m is filled if and only if it does not have D_10 as a quotient.")
    print("This condition forces the group G to be a direct product G = Q x C_2, where Q is a group of order 5^m.")
    print("For G to be nonabelian, Q must be a nonabelian group.")
    print("Groups of order p^m (p prime) are only nonabelian if m >= 3 (for p odd).")
    print(f"Thus, for q = 5, the nonabelian filled groups are those of the form Q x C_2, where Q is any nonabelian group of order 5^m. This requires m >= 3.")

    # Case 3: q >= 7
    print("\nCase 3: q >= 7 is a prime")
    print("The order of the group is 2 * q^m, where q >= 7.")
    print("The prime factors of the order of G are 2 and q.")
    print("As noted, the prime factors of the orders of all forbidden groups are from {2, 3, 5}.")
    print("Since q is not 3 or 5, the order of any forbidden group cannot divide 2 * q^m (unless the forbidden group's order is a power of 2, like 8, which does not divide 2*q^m for odd q).")
    print("This means a group of order 2*q^m (for q>=7) can never have a forbidden group as a quotient.")
    print("Therefore, for q >= 7, ALL nonabelian groups of order 2 * q^m are filled groups for any m >= 1.")

    print("\n--- Summary of Results ---")
    print("The nonabelian filled groups of order 2*q^m are as follows:")
    print("1. If q = 3: There are no such groups.")
    print("2. If q = 5: The groups of the form Q x C_2, where Q is a nonabelian group of order 5^m. This requires m >= 3.")
    print("3. If q >= 7: All nonabelian groups of order 2*q^m are filled.")


if __name__ == '__main__':
    solve_filled_groups_problem()
    print("\n<<<")
    print("q = 3: None.")
    print("q = 5: Groups of the form Q x C_2 where Q is a nonabelian group of order 5^m, for m >= 3.")
    print("q >= 7: All nonabelian groups of order 2*q^m.")
    print(">>>")
