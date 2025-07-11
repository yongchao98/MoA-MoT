def solve_group_theory_problem():
    """
    This function finds the minimum value of y (number of Sylow 5-subgroups)
    that guarantees a group G is nonsolvable, given n_3 (number of Sylow
    3-subgroups) is at most 9.

    The solution is based on established theorems in finite group theory.
    """
    
    # According to Sylow's third theorem, n_p must be congruent to 1 mod p.
    # We search for the smallest y such that y % 5 == 1.
    y = 1
    while True:
        # Check if y is a possible number of Sylow 5-subgroups.
        if y % 5 == 1:
            # Check if this y guarantees nonsolvability.
            # Case y = 1 (n_5 = 1):
            # A group with a normal Sylow subgroup is not necessarily solvable.
            # However, we can easily find a solvable group that satisfies the conditions.
            # Example: G = C_15 (cyclic group of order 15).
            # n_3(C_15) = 1 (<= 9).
            # n_5(C_15) = 1.
            # C_15 is abelian, hence solvable.
            # So, y=1 does not guarantee nonsolvability.
            if y == 1:
                y += 1
                continue

            # Case y > 1. The next possible value is y=6.
            # Case y = 6 (n_5 = 6):
            # It's a known theorem in finite group theory that if a group G has exactly 6 Sylow 5-subgroups,
            # then G has a quotient group G/N which is isomorphic to A_5 or S_5.
            # Both A_5 and S_5 are nonsolvable groups.
            # If a group has a nonsolvable quotient, it must be nonsolvable itself.
            # Therefore, if n_5 = 6, the group G must be nonsolvable.
            # This conclusion holds irrespective of the number of Sylow 3-subgroups.
            # The condition n_3 <= 9 is extra information that does not alter the conclusion for y=6.
            
            # Since y=1 does not work, and y=6 is the next possible value that does work,
            # the minimum value of y is 6.
            
            final_y = y

            print(f"Let y be the number of Sylow 5-subgroups, n_5.")
            print(f"Let n_3 be the number of Sylow 3-subgroups.")
            print(f"We are given n_3 <= 9 and n_5 = y.")
            print(f"We need to find the minimum y that forces G to be nonsolvable.")
            print(f"By Sylow's Theorems, y must satisfy y = 1 (mod 5).")
            print(f"Possible values for y are 1, 6, 11, 16, ...")
            print("")
            print(f"Step 1: Test y = 1.")
            print(f"If y = 1, consider the cyclic group G = C_15. |G|=15.")
            print(f"For C_15, n_3 = 1 (<= 9) and n_5 = 1.")
            print(f"C_15 is abelian, and therefore solvable. So y = 1 does not guarantee nonsolvability.")
            print("")
            print(f"Step 2: Test y = 6.")
            print(f"If y = n_5 = 6, G acts on its 6 Sylow 5-subgroups by conjugation.")
            print(f"This induces a homomorphism phi: G -> S_6.")
            print(f"The quotient group G/ker(phi) is isomorphic to a transitive subgroup of S_6 containing an element of order 5.")
            print(f"By classification of such subgroups, G/ker(phi) is isomorphic to A_5 or S_5.")
            print(f"Since A_5 and S_5 are nonsolvable, G must be nonsolvable.")
            print(f"This conclusion holds regardless of the value of n_3.")
            print("")
            print(f"Step 3: Conclusion.")
            print(f"y=1 is not sufficient. y=6 is sufficient.")
            print(f"The minimum value for y is {final_y}.")
            return final_y

        y += 1

solve_group_theory_problem()