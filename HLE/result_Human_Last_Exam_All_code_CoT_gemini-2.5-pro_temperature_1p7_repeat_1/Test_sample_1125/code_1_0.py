def solve_hat_puzzle():
    """
    Calculates the maximum number of people guaranteed to know their hat number.
    The strategy is a multi-round tournament-style elimination.
    """
    
    total_people = 12
    unsolved_people = total_people

    # --- Strategy Explanation ---
    print("The team devises a 'tournament' strategy to solve the puzzle.")
    print(f"Initially, there are {unsolved_people} people who don't know their number.")
    print("-" * 30)

    # Round 1
    # 6 pairs are formed from 12 people. 6 revelations are made.
    round1_solved = unsolved_people // 2
    unsolved_people -= round1_solved
    print(f"Round 1: {total_people} people form {round1_solved} pairs.")
    print(f"           One number is revealed per pair.")
    print(f"           People who know their number: {round1_solved}")
    print(f"           Remaining unsolved people: {unsolved_people}")
    print("-" * 30)

    # Round 2
    # 3 pairs are formed from the 6 unsolved people. 3 revelations are made.
    round2_solved = unsolved_people // 2
    unsolved_people -= round2_solved
    print(f"Round 2: The {round2_solved*2} unsolved people form {round2_solved} pairs.")
    print(f"           One number is revealed per pair.")
    print(f"           Additional people who know their number: {round2_solved}")
    print(f"           Remaining unsolved people: {unsolved_people}")
    print("-" * 30)

    # Round 3
    # 1 pair is formed from the 3 unsolved people. 1 revelation is made.
    round3_solved = unsolved_people // 2
    unsolved_people -= round3_solved
    print(f"Round 3: From the {round3_solved*2+1} unsolved people, {round3_solved} pair is formed.")
    print(f"           One number is revealed from the pair.")
    print(f"           Additional people who know their number: {round3_solved}")
    print(f"           Remaining unsolved people: {unsolved_people}")
    print("-" * 30)

    # Round 4
    # The final pair is formed. 1 revelation is made.
    round4_solved = unsolved_people // 2
    unsolved_people -= round4_solved
    print(f"Round 4: The final {round4_solved*2} unsolved people form a pair.")
    print(f"           One number is revealed from the pair.")
    print(f"           Additional people who know their number: {round4_solved}")
    print(f"           Remaining unsolved people: {unsolved_people}")
    print("-" * 30)

    # Final Deduction
    # The last person deduces their number.
    final_deduction_solved = unsolved_people
    print("Final Step: The last person's number was never revealed directly.")
    print("However, they know the other 11 unique numbers and can deduce their own by elimination.")
    print(f"People solved by deduction: {final_deduction_solved}")
    print("-" * 30)

    # Final Calculation
    total_solved = round1_solved + round2_solved + round3_solved + round4_solved + final_deduction_solved
    
    print("The total number (N) of people guaranteed to determine their number is the sum of those solved in each step:")
    print(f"{round1_solved} + {round2_solved} + {round3_solved} + {round4_solved} + {final_deduction_solved} = {total_solved}")

solve_hat_puzzle()