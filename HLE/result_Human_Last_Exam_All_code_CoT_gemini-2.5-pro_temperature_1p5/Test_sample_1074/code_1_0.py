import sys

def solve_group_theory_problem():
    """
    Solves for the minimum value of y based on the given group theory conditions.
    This function programmatically explains the logical steps to reach the conclusion.
    """
    
    # Step 1: Define the problem
    print("Goal: Find the minimum integer y such that if G is a finite group with:")
    print("  1. The number of Sylow 3-subgroups (n_3) is at most 9,")
    print("  2. The number of Sylow 5-subgroups (n_5) is y,")
    print("then G must be nonsolvable.")
    print("-" * 40)

    # Step 2: Analyze constraints from Sylow's Theorems
    print("Sylow's third theorem states that n_p must be congruent to 1 modulo p.")
    
    # Constraints on n_3
    possible_n3 = [n for n in range(1, 10) if n % 3 == 1]
    print(f"From n_3 <= 9 and n_3 % 3 == 1, the possible values for n_3 are: {possible_n3}.")

    # Constraints on n_5 = y
    possible_y = [n for n in range(1, 20) if n % 5 == 1]
    print(f"From y % 5 == 1, the possible values for y are {possible_y}...")
    print("-" * 40)

    # Step 3: Test the smallest possible values of y
    print("We test the possible values for y, starting from the smallest.")
    
    # Test y = 1
    y_test = 1
    print(f"\n--- Testing y = {y_test} ---")
    print(f"If y = {y_test}, then n_5 = 1. A group with n_p = 1 has a normal Sylow p-subgroup.")
    print("Let's check if there is a SOLVABLE group G where n_3 <= 9 and n_5 = 1.")
    print("Consider the group G = Z_15 (the cyclic group of order 15).")
    print("For G = Z_15:")
    print(" - n_3 = 1, which satisfies n_3 <= 9.")
    print(" - n_5 = 1.")
    print(" - Z_15 is abelian, and therefore solvable.")
    print(f"Conclusion: Since a solvable group exists for y = {y_test}, this value is not the answer.")
    
    # Test y = 6
    y_test = 6
    print(f"\n--- Testing y = {y_test} ---")
    print(f"This is the next smallest possible value for y since {y_test} % 5 == 1.")
    print("Let's prove that if any group G has n_5 = 6, it must be nonsolvable.")
    print("1. Let G act on its set of 6 Sylow 5-subgroups by conjugation.")
    print("2. This action gives a homomorphism phi: G -> S_6 (symmetric group on 6 elements).")
    print("3. Let K = ker(phi). Then G/K is isomorphic to a transitive subgroup of S_6.")
    print("4. The order of G/K must be divisible by n_5 = 6.")
    print("5. A Sylow 5-subgroup P cannot be in K (otherwise n_5 would be 1).")
    print("   Thus, the order of G/K must also be divisible by 5.")
    print("6. Therefore, the order of G/K must be divisible by lcm(5, 6) = 30.")

    print("\n   The transitive subgroups of S_6 with order divisible by 30 are A_5, S_5, A_6, and S_6.")
    print("   - |A_5| = 60")
    print("   - |S_5| = 120")
    print("   - |A_6| = 360")
    print("   - |S_6| = 720")
    print("\n7. All of these possible groups for G/K are nonsolvable (as they contain a non-abelian simple group, A_5 or A_6).")
    print("8. If a group G has a nonsolvable quotient G/K, G itself must be nonsolvable.")
    
    print("\nConclusion: The condition n_5 = 6 ALONE forces G to be nonsolvable.")
    print("This means that if (n_3 <= 9 AND n_5 = 6), G must be nonsolvable.")
    print("-" * 40)

    # Step 4: Final Answer
    final_y = 6
    print("We have shown:")
    print("- y = 1 does not guarantee nonsolvability.")
    print("- y = 6, the next possible value, DOES guarantee nonsolvability.")
    print("\nTherefore, the minimum value of y is 6.")
    print(f"\nThe final equation can be stated as: y = {final_y}")

solve_group_theory_problem()