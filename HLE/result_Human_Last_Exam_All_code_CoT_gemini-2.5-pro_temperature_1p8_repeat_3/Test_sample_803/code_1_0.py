def solve_group_theory_question():
    """
    This function provides the solution to the group theory question about nonabelian filled groups.
    It prints the theoretical characterization of these groups along with examples.
    """

    print("The nonabelian filled groups of order 2q^m for an odd prime q and a natural number m are characterized as follows:\n")
    print("A group G is a nonabelian filled group of order 2q^m if and only if it can be expressed as a semidirect product:")
    print("\n    G = Q \u22ca C\u2082\n")
    print("where:\n")
    print(f"1. Q is a group of order q^m.")
    print(f"2. C\u2082 is the cyclic group of order 2.")
    print(f"3. The action of C\u2082 on Q is given by an automorphism \u03B1 of Q. This automorphism \u03B1 must satisfy the following conditions:")
    print(f"   a) \u03B1 is not the trivial automorphism (which ensures G is nonabelian).")
    print(f"   b) The action of \u03B1 on the abelianization of Q (denoted Q/Q') is the inversion map.")
    print(f"\nThe second condition on the automorphism \u03B1 means that for every element g in Q, the following relation holds:")
    print(f"\n    \u03B1(g) * g \u2208 Q' \n")
    print(f"where Q' is the commutator subgroup of Q. This condition ensures the group G is 'filled'.\n")

    print("--- Examples ---\n")

    # Example 1: Order 6 (q=3, m=1)
    print("Example for order 6 (q=3, m=1):")
    print("Let Q = C\u2083, the cyclic group of order 3. Q is abelian, so its commutator Q' is the trivial group {e}.")
    print("The automorphism \u03B1 must act as inversion on Q/Q' = Q = C\u2083. The inversion map \u03B1(x) = x\u207b\u00b9 is a valid automorphism.")
    print("The resulting group G = C\u2083 \u22ca C\u2082 is isomorphic to the dihedral group D\u2086 (also the symmetric group S\u2083). This is the only nonabelian group of order 6, and it is filled.\n")

    # Example 2: Order 10 (q=5, m=1)
    print("Example for order 10 (q=5, m=1):")
    print("Let Q = C\u2085. Q' = {e}. \u03B1 is the inversion map on C\u2085.")
    print("The resulting group G = C\u2085 \u22ca C\u2082 is the dihedral group D\u2081\u2080. This is the only nonabelian group of order 10, and it is filled.\n")

    # Example 3: Order 18 (q=3, m=2)
    print("Example for order 18 (q=3, m=2):")
    print("There are two groups Q of order 9 = 3\u00b2, both are abelian: C\u2089 and C\u2083 x C\u2083.")
    print("1. If Q = C\u2089, its commutator Q'={e}. Let \u03B1 be the inversion map. The group G = C\u2089 \u22ca C\u2082 is the dihedral group D\u2081\u2088.")
    print("2. If Q = C\u2083 x C\u2083, its commutator Q'={e}. Let \u03B1 be the inversion map. The group G = (C\u2083 x C\u2083) \u22ca C\u2082.")
    print("Both of these nonabelian groups of order 18 are filled.\n")
    
    # Example 4: Order 54 (q=3, m=3, Q nonabelian)
    print("Example for order 54 (q=3, m=3) with a nonabelian Q:")
    print("Let Q be the nonabelian group of order 27 and exponent 3.")
    print("Its commutator subgroup Q' is isomorphic to C\u2083, and its abelianization Q/Q' is isomorphic to C\u2083 x C\u2083.")
    print("It can be shown that there exists an automorphism \u03B1 of Q which acts as inversion on Q/Q'.")
    print("The resulting semidirect product G = Q \u22ca C\u2082 is a nonabelian filled group of order 54.")

if __name__ == '__main__':
    solve_group_theory_question()