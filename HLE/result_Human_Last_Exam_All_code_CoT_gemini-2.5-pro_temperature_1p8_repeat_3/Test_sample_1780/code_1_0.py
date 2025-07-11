import math

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters L, N, S, and W
    based on the described connection rule.
    """
    # Based on our analysis, the problem reduces to a permutation calculation.
    # The letters that can appear in the first three positions are L, N, and S.
    # 'W' must be the last letter.
    # Therefore, we need to find the number of ways to arrange {L, N, S}.
    num_letters_to_permute = 3

    # The number of permutations of n distinct items is n! (n factorial).
    result = math.factorial(num_letters_to_permute)

    print(f"The problem reduces to finding the number of permutations of the letters L, N, and S, which must be followed by W.")
    print(f"The number of ways to arrange {num_letters_to_permute} distinct letters is {num_letters_to_permute}! (factorial).")
    print(f"The final equation is: 3 * 2 * 1 = {result}")
    print(f"\nThus, the total number of valid arrangements is {result}.")

solve_letter_arrangement()