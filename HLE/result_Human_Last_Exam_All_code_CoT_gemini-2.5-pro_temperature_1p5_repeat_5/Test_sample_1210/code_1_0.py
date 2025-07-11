def solve_mancala_difference():
    """
    Calculates which score differences are possible in a Mancala game
    with a specific total number of stones.
    """
    # Step 1: Calculate the total number of stones from the initial state.
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21

    total_stones = sum(p1_pits) + p1_store + sum(p2_pits) + p2_store

    print(f"The total number of stones on the board is {total_stones}.")
    print("At the end of the game, let the winner's score be S_winner and the loser's score be S_loser.")
    print(f"We know that S_winner + S_loser = {total_stones}.")
    print("The score difference D is defined as S_winner - S_loser.")
    print("By solving for the winner's score, we get the formula: S_winner = (total_stones + D) / 2.")
    print("For S_winner to be a whole number, (total_stones + D) must be an even number.")
    print("-" * 20)

    # Step 2: Define the answer choices and check each one.
    answer_choices = {
        'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5
    }
    unobtainable_choices = []

    for label, diff in answer_choices.items():
        numerator = total_stones + diff
        # The winner's score must be an integer.
        is_possible = (numerator % 2 == 0)
        
        winner_score = numerator / 2.0
        
        print(f"For answer choice {label}, the score difference (D) is {diff}.")
        # Final code outputs each number in the final equation as requested.
        print(f"The winner's score would be: ({total_stones} + {diff}) / 2 = {winner_score}")

        if is_possible:
            print("This is a possible score difference because the winner's score is a whole number.")
        else:
            print("This is NOT a possible score difference because the winner's score is not a whole number.")
            unobtainable_choices.append(label)
        print()

    # Step 3: Conclude based on the number of unobtainable differences.
    print("-" * 20)
    if len(unobtainable_choices) > 1:
        print(f"The score differences for choices {', '.join(unobtainable_choices)} are unobtainable.")
        print("Therefore, more than one of the listed score differences is unobtainable.")
    elif len(unobtainable_choices) == 1:
        print(f"The score difference for choice {unobtainable_choices[0]} is the only one that is unobtainable.")
    else:
        print("All listed score differences are obtainable.")

solve_mancala_difference()
<<<G>>>