def find_minimum_y():
    """
    This script determines the minimum value of y (number of Sylow 5-subgroups)
    that forces a group G to be nonsolvable, given n_3 (number of Sylow 3-subgroups) <= 9.
    """
    
    # Step 1: Analyze the constraints on n_3 and n_5 from Sylow's Theorems.
    print("Step 1: Analyze the constraints from Sylow's Theorems.")
    print("Sylow's Third Theorem states that for a prime p, the number of Sylow p-subgroups (n_p) must be congruent to 1 modulo p (n_p % p == 1).")
    
    # For n_3:
    n3_limit = 9
    print(f"\nFor p=3, n_3 must be congruent to 1 mod 3. The condition is n_3 <= {n3_limit}.")
    print("Possible values for n_3 are 1, 4, 7, 10, ...")
    print(f"Given n_3 <= {n3_limit}, the allowed values for n_3 are {{1, 4, 7}}.")

    # For n_5 = y:
    print("\nFor p=5, n_5 (or y) must be congruent to 1 mod 5.")
    print("Possible values for y are {1, 6, 11, 16, 21, ...}.")
    print("-" * 30)

    # Step 2: Test the smallest possible value for y, which is 1.
    print("Step 2: Test the smallest possible value, y = 1.")
    y_test_1 = 1
    print(f"If y = {y_test_1}, the conditions are n_3 <= {n3_limit} and n_5 = {y_test_1}.")
    print("To check if this forces nonsolvability, we see if a SOLVABLE group can satisfy these conditions.")
    print("Let's consider the group G = Z_5 x S_3 (the direct product of the cyclic group of order 5 and the symmetric group of degree 3).")
    print("G is a direct product of two solvable groups (Z_5 is abelian, S_3 is solvable), so G is solvable.")
    print("Let's check its n_3 and n_5:")
    print("n_5(G) = n_5(Z_5) * n_5(S_3) = 1 * 1 = 1. This satisfies n_5 = 1.")
    print("n_3(G) = n_3(Z_5) * n_3(S_3) = 1 * 1 = 1. This satisfies n_3 <= 9.")
    print("\nConclusion for y=1: A solvable group exists that meets the criteria. Therefore, y=1 does not force nonsolvability.")
    print("-" * 30)

    # Step 3: Test the next possible value for y, which is 6.
    print("Step 3: Test the next possible value, y = 6.")
    y_test_2 = 6
    print(f"If y = {y_test_2}, the conditions are n_3 <= {n3_limit} and n_5 = {y_test_2}.")
    print("This forces nonsolvability if NO solvable group with n_5 = 6 exists.")
    
    print("\nWe use a known theorem from finite group theory:")
    print("Theorem: If a group G has n_p = p + 1 Sylow p-subgroups, then G/O_p(G) is isomorphic to A_{p+1} or S_{p+1}.")
    print("(O_p(G) is the largest normal p-subgroup of G).")
    
    p = 5
    n_p_val = p + 1
    print(f"\nIn our case, p = {p} and y = n_5 = {y_test_2}, which matches the theorem's condition n_p = p + 1 = {n_p_val}.")
    print(f"The theorem implies that for any group G with n_5 = 6, the quotient group G/O_5(G) must be isomorphic to A_6 or S_6.")
    
    print("\nNow, assume G is a solvable group.")
    print("A fundamental property of solvable groups is that any quotient group of a solvable group is also solvable.")
    print("This means G/O_5(G) must be solvable.")
    print("However, A_6 is a non-abelian simple group, making it nonsolvable. S_6 contains A_6 as a normal subgroup, so it is also nonsolvable.")
    print("This is a contradiction. The assumption that G is solvable must be false.")
    
    print("\nConclusion for y=6: Any group with n_5 = 6 must be nonsolvable. This is true regardless of the value of n_3.")
    print("-" * 30)

    # Step 4: Final determination of the minimum value.
    print("Step 4: Determine the minimum value of y.")
    print(f"We found that y = {y_test_1} does not force nonsolvability.")
    print(f"We found that y = {y_test_2} *does* force nonsolvability.")
    print(f"Since {y_test_2} is the next possible value for y after {y_test_1}, it is the minimum value required.")
    
    min_y = 6
    print("\n" + "="*40)
    print("THE FINAL ANSWER")
    print("="*40)
    print(f"The minimum value of y is {min_y}.")
    print("The final implication with this minimum value is:")
    print(f"If the number of Sylow 3-subgroups is at most {n3_limit} and the number of Sylow 5-subgroups is {min_y}, then the group G must be nonsolvable.")

find_minimum_y()