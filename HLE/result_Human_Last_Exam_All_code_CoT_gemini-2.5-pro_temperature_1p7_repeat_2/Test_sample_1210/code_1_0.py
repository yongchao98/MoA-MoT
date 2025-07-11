def solve_mancala_difference():
    """
    Calculates the possible score differences in a Mancala game
    based on the parity of the total number of stones.
    """
    # Step 1: Define the game state and calculate the total number of stones.
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21

    sum_p1_pits = sum(p1_pits)
    sum_p2_pits = sum(p2_pits)
    
    total_stones = sum_p1_pits + p1_store + sum_p2_pits + p2_store

    print(f"The total number of stones (T) is calculated as: {sum_p1_pits} + {p1_store} + {sum_p2_pits} + {p2_store} = {total_stones}")

    # Step 2: Explain the logic based on parity.
    print("\nAt the end of any Mancala game, let the winner's score be W and the loser's score be L.")
    print(f"The sum of scores must equal the total stones: W + L = T = {total_stones}")
    print("The score difference is D: W - L = D")
    
    # Step 3: Show the final equation with numbers.
    print("\nBy adding these two equations, we get: 2 * W = T + D")
    print(f"For this specific game, the equation is: 2 * W = {total_stones} + D")
    print("\nFor W to be a whole number, the right side of the equation must be an even number.")
    print(f"Since the total number of stones ({total_stones}) is even, the score difference (D) must also be even.")

    # Step 4: Check the given answer choices.
    answer_choices = { 'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5 }
    untainable_differences = []

    print("\nAnalyzing the score differences from the answer choices:")
    for key, diff in answer_choices.items():
        if diff % 2 != 0:
            untainable_differences.append(diff)
            print(f"- Choice {key}: A difference of {diff} is an ODD number, which is mathematically unobtainable.")
        else:
            print(f"- Choice {key}: A difference of {diff} is an EVEN number, which is potentially obtainable.")
    
    # Step 5: Formulate the final conclusion.
    print(f"\nThe unobtainable score differences are {untainable_differences}.")
    if len(untainable_differences) > 1:
        print("Conclusion: More than one of the listed score differences is unobtainable.")
    elif len(untainable_differences) == 1:
        untainable_key = [k for k, v in answer_choices.items() if v == untainable_differences[0]][0]
        print(f"Conclusion: The only unobtainable score difference is {untainable_differences[0]} (Choice {untainable_key}).")
    else:
        print("Conclusion: All listed score differences are potentially obtainable.")

solve_mancala_difference()