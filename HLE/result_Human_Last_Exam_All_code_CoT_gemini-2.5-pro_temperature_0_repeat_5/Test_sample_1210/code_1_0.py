import math

def solve_mancala_difference():
    """
    Analyzes a Mancala game state to determine which score differences are impossible.
    """
    # Initial game state
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21

    # Step 1: Calculate the total number of stones in the game.
    total_stones = sum(p1_pits) + p1_store + sum(p2_pits) + p2_store

    print("Step 1: Calculate the total number of stones.")
    print("The total number of stones is the sum of all stones in all pits and stores.")
    # The prompt asks to output each number in the final equation.
    # Here is the equation for the total number of stones.
    print(f"Total = ({' + '.join(map(str, p1_pits))}) + {p1_store} + ({' + '.join(map(str, p2_pits))}) + {p2_store}")
    print(f"Total = {sum(p1_pits)} + {p1_store} + {sum(p2_pits)} + {p2_store}")
    print(f"Total = {sum(p1_pits) + p1_store} + {sum(p2_pits) + p2_store}")
    print(f"Total = {total_stones}")
    print("-" * 40)

    # Step 2: Analyze the properties of the final score difference.
    print("Step 2: Analyze the relationship between scores and the total.")
    print("Let S1 be the final score of Player 1 and S2 be the final score of Player 2.")
    print("At the end of the game, the sum of the scores must equal the total number of stones.")
    print(f"Equation 1: S1 + S2 = {total_stones}")
    print()
    print("The score difference, D, is the winner's score minus the loser's score.")
    print("Equation 2: |S1 - S2| = D")
    print()
    print("Let's solve these equations for S1 (assuming S1 is the winner's score):")
    print(f"From (S1 + S2 = {total_stones}) and (S1 - S2 = D), we can add them together:")
    print(f"(S1 + S2) + (S1 - S2) = {total_stones} + D")
    print(f"2 * S1 = {total_stones} + D")
    # Here is the final equation with each number.
    print(f"S1 = ({total_stones} + D) / 2  =>  S1 = {int(total_stones / 2)} + D / 2")
    print()
    print(f"Since S1 (a player's score) must be a whole number, and {int(total_stones / 2)} is a whole number,")
    print("then D/2 must also be a whole number. This is only true if D is an EVEN number.")
    print("Therefore, any odd score difference is impossible.")
    print("-" * 40)

    # Step 3: Evaluate the answer choices.
    print("Step 3: Evaluate the given answer choices.")
    answer_choices = { "A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5 }
    impossible_choices = []

    for choice, value in answer_choices.items():
        if value % 2 != 0:
            print(f" - Score difference {value} (Choice {choice}) is ODD, so it is impossible.")
            impossible_choices.append(choice)

    print()
    if len(impossible_choices) > 1:
        print(f"We have found {len(impossible_choices)} impossible score differences in the list ({', '.join(impossible_choices)}).")
        print("This means that 'More than one of the listed score differences is unobtainable' is the correct statement.")
    else:
        # This case will not be reached based on the logic.
        print("Only one impossible score difference was found.")

solve_mancala_difference()