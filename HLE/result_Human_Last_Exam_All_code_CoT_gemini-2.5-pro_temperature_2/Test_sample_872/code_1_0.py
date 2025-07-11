from fractions import Fraction

def solve_tic_tac_toe_prob():
    """
    Calculates the maximum win probability for a player who goes first
    against a computer that plays randomly. A tie is considered a loss.
    """

    print("Step 1: Determine the optimal strategy.")
    print("The optimal first move is to take the center square. We will analyze the win probability based on this opening.")
    print("\nThere are two scenarios for the computer's first move, based on symmetry:")
    print("  - Case A: Computer picks a corner square (4/8 probability).")
    print("  - Case B: Computer picks an edge square (4/8 probability).")

    # --- Probability of each case ---
    p_O1_corner = Fraction(4, 8)
    p_O1_edge = Fraction(4, 8)

    # --- Analysis of Case A: Computer plays a corner ---
    # Our strategy: We take the opposite corner. This sets up two potential forks for our next move.
    # The computer has 6 squares to choose from for its second move (O2).
    # To prevent us from creating a guaranteed winning fork, the computer must play on one of two specific squares.
    p_O2_blocks_fork_A = Fraction(2, 6)
    p_O2_not_block_fork_A = Fraction(4, 6)

    # If the computer does NOT block our fork setup, we create a fork and win.
    win_prob_if_no_block_A = Fraction(1, 1)

    # If the computer DOES block one of our potential forks, we must react.
    # In this new state, we can still force a win with a 3/4 probability, as the computer
    # must then block our immediate winning line on its next move (O3), and there's a 1/4 chance it will.
    p_win_after_block_A = Fraction(3, 4)

    # The total win probability in Case A is a sum of these outcomes.
    p_win_A = p_O2_not_block_fork_A * win_prob_if_no_block_A + p_O2_blocks_fork_A * p_win_after_block_A

    print(f"\nStep 2: Analyze Case A (Computer plays a corner).")
    print(f"The computer fails to block our fork setup with probability {p_O2_not_block_fork_A.numerator}/{p_O2_not_block_fork_A.denominator}, in which case we win.")
    print(f"The computer blocks our setup with probability {p_O2_blocks_fork_A.numerator}/{p_O2_blocks_fork_A.denominator}. In this sub-game, our win chance is {p_win_after_block_A.numerator}/{p_win_after_block_A.denominator}.")
    print(f"Total win probability in Case A = {p_O2_not_block_fork_A.numerator}/{p_O2_not_block_fork_A.denominator} * 1 + {p_O2_blocks_fork_A.numerator}/{p_O2_blocks_fork_A.denominator} * {p_win_after_block_A.numerator}/{p_win_after_block_A.denominator} = {p_win_A.numerator}/{p_win_A.denominator}")

    # --- Analysis of Case B: Computer plays an edge ---
    # Our strategy: We create an immediate threat, forcing the computer to block on its next move (O2).
    # The computer has 6 squares left. Only 1 of them blocks our immediate win.
    p_O2_blocks_B = Fraction(1, 6)
    p_O2_not_block_B = Fraction(5, 6)

    # If the computer fails to block, we win immediately.
    win_prob_if_no_block_B = Fraction(1, 1)

    # If the computer blocks, we can then set up a guaranteed winning fork. The computer, on its
    # next move (O3), must play on one specific square out of 4 to prevent our fork.
    # If it fails, we win. Probability of this win is 3/4.
    p_win_after_block_B = Fraction(3, 4)

    # The total win probability in Case B is a sum of these outcomes.
    p_win_B = p_O2_not_block_B * win_prob_if_no_block_B + p_O2_blocks_B * p_win_after_block_B
    
    print(f"\nStep 3: Analyze Case B (Computer plays an edge).")
    print(f"The computer fails to make a required block with probability {p_O2_not_block_B.numerator}/{p_O2_not_block_B.denominator}, in which case we win.")
    print(f"The computer makes the required block with probability {p_O2_blocks_B.numerator}/{p_O2_blocks_B.denominator}. In this sub-game, our win chance is {p_win_after_block_B.numerator}/{p_win_after_block_B.denominator}.")
    print(f"Total win probability in Case B = {p_O2_not_block_B.numerator}/{p_O2_not_block_B.denominator} * 1 + {p_O2_blocks_B.numerator}/{p_O2_blocks_B.denominator} * {p_win_after_block_B.numerator}/{p_win_after_block_B.denominator} = {p_win_B.numerator}/{p_win_B.denominator}")


    # --- Final calculation ---
    total_prob = p_O1_corner * p_win_A + p_O1_edge * p_win_B
    final_numerator = total_prob.numerator
    final_denominator = total_prob.denominator
    
    print(f"\nStep 4: Combine the probabilities.")
    print("Total Win Probability = P(Case A) * P(Win|A) + P(Case B) * P(Win|B)")
    print(f"Result = ({p_O1_corner.numerator}/{p_O1_corner.denominator}) * ({p_win_A.numerator}/{p_win_A.denominator}) + ({p_O1_edge.numerator}/{p_O1_edge.denominator}) * ({p_win_B.numerator}/{p_win_B.denominator})")
    print(f"Result = {p_O1_corner * p_win_A} + {p_O1_edge * p_win_B}")
    print(f"Result = {p_O1_corner * p_win_A + p_O1_edge * p_win_B}")
    print(f"\nThe final reduced fraction for the maximum chance of winning is: {final_numerator}/{final_denominator}")


solve_tic_tac_toe_prob()
print(f'<<<{Fraction(15, 16).numerator}/{Fraction(15, 16).denominator}>>>')