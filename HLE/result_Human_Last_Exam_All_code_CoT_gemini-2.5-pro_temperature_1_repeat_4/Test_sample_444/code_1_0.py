import random

def strategy_A_simulation(true_sequence):
    """
    Simulates Alice's randomized strategy for Case (A).
    This strategy is based on guessing the modulo of the sum of all numbers.
    It demonstrates why it's hard to achieve the desired success rate.
    """
    print(f"The true sequence is S = {true_sequence}...")

    # Alice's plan:
    # 1. She will leave box 1 closed and open all other boxes.
    # 2. She will calculate the sum of the numbers in the opened boxes.
    # 3. She will pick a random 'j' from {0, ..., 9} and assume the sum of ALL numbers is 'j' modulo 10.
    # 4. Based on this assumption, she will calculate her one guess for the number in box 1.

    closed_box_index = 1
    opened_boxes_values = true_sequence[1:]
    true_value_in_box_1 = true_sequence[0]

    # Alice computes the sum of the numbers she sees. This is well-defined as the sequence is eventually zero.
    sum_of_opened_values = sum(opened_boxes_values)

    # Alice makes a random choice for what she thinks the total sum will be modulo 10.
    j = random.randint(0, 9)
    print(f"Alice randomly assumes the total sum modulo 10 will be j = {j}.")

    # Alice calculates her guess 'g'. Her guess must be the number that makes her assumption true.
    # She guesses the smallest non-negative integer 'g' that satisfies:
    # (g + sum_of_opened_values) % 10 == j
    sum_mod_10 = sum_of_opened_values % 10
    g = (j - sum_mod_10 + 10) % 10

    print(f"Alice opens all boxes except box 1 and sees the sum of these is {sum_of_opened_values}.")
    print("Her guess 'g' for the number in box 1 is calculated from the equation:")
    print(f"(g + {sum_of_opened_values}) mod 10 = {j}")
    print(f"The calculated guess is g = {g}.")
    
    print("-" * 20)
    print(f"Alice's Final Guess for Box 1: {g}")
    print(f"The True Value in Box 1: {true_value_in_box_1}")

    if g == true_value_in_box_1:
        print("Result: SUCCESS")
    else:
        print("Result: FAILURE")
    
    # Analysis of the success probability
    true_total_sum = sum(true_sequence)
    true_sum_mod_10 = true_total_sum % 10
    print("\n--- Analysis ---")
    print("This strategy succeeds if Alice's random choice 'j' matches the true value of the total sum modulo 10.")
    print(f"The true total sum is {true_total_sum}, so the true value modulo 10 is {true_sum_mod_10}.")
    print(f"Since Alice picks 'j' randomly from 10 possibilities, her chance of success is 1/10 for any sequence.")
    print("A 1/10 chance of success is not >= 9/10.")

# Let's run the simulation. The adversary can choose any eventually-zero sequence.
adversary_sequence = [1, 8, 4, 0, 2, 0, 0, 0]
strategy_A_simulation(adversary_sequence)