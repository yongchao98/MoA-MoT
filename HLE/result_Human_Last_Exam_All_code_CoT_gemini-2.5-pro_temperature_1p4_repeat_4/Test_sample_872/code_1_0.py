from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum win probability for a player against a random-move computer
    in Tic Tac Toe, where the player goes first.
    """

    # Step 1 & 2: The player's optimal first move is the center.
    # The computer has 8 squares left. We analyze its two possible move types.
    p_c1_corner = Fraction(4, 8)
    p_c1_edge = Fraction(4, 8)

    print("Player starts by placing 'X' in the center square.")
    print(f"The computer can either pick a corner (Prob={p_c1_corner}) or an edge (Prob={p_c1_edge}).\n")

    # After P1, there are 8 squares. After C1, there are 7. After P2, there are 6.
    # P2 will always create a line of 'X-X-.'
    # C2 has 6 squares to choose from. Only 1 of them blocks the player's immediate win.
    p_c2_blocks = Fraction(1, 6)
    p_c2_fails_to_block = Fraction(5, 6)
    
    # --- Case A: Computer's first move is a corner ---
    print("--- Analysis for Case A: Computer plays a corner first ---")
    # If C2 blocks P2's threat, we need to find the win probability from this new state.
    # The board state is: O(corner), X(center), X(adj edge), O(block).
    # P3 can create another threat that C3 will fail to block with probability 2/3.
    # If C3 also blocks (prob 1/3), the game results in a draw.
    win_prob_if_c2_blocks_in_corner_case = Fraction(2, 3)
    
    p_win_given_c1_corner = p_c2_fails_to_block * 1 + p_c2_blocks * win_prob_if_c2_blocks_in_corner_case

    print("Player's second move creates an immediate threat.")
    print(f"The computer fails to block this threat with probability {p_c2_fails_to_block}, and the player wins.")
    print(f"The computer blocks this threat with probability {p_c2_blocks}.")
    print(f"If the computer blocks, the player can make another move that wins with probability {win_prob_if_c2_blocks_in_corner_case}.")
    print(f"So, the total win probability if the computer first plays a corner is:")
    print(f"P(Win|C1=corner) = {p_c2_fails_to_block}*1 + {p_c2_blocks}*{win_prob_if_c2_blocks_in_corner_case} = {p_win_given_c1_corner}\n")

    # --- Case B: Computer's first move is an edge ---
    print("--- Analysis for Case B: Computer plays an edge first ---")
    # If C2 blocks P2's threat, the resulting game state forces the player to be on the defensive
    # and the game leads to a draw if both sides play correctly. A draw is a loss for the player.
    win_prob_if_c2_blocks_in_edge_case = Fraction(0, 1)

    p_win_given_c1_edge = p_c2_fails_to_block * 1 + p_c2_blocks * win_prob_if_c2_blocks_in_edge_case

    print("Player's second move creates an immediate threat.")
    print(f"The computer fails to block this threat with probability {p_c2_fails_to_block}, and the player wins.")
    print(f"The computer blocks this threat with probability {p_c2_blocks}.")
    print(f"If the computer blocks, the game leads to a draw (win probability = {win_prob_if_c2_blocks_in_edge_case}).")
    print(f"So, the total win probability if the computer first plays an edge is:")
    print(f"P(Win|C1=edge) = {p_c2_fails_to_block}*1 + {p_c2_blocks}*{win_prob_if_c2_blocks_in_edge_case} = {p_win_given_c1_edge}\n")

    # --- Final Calculation ---
    total_win_prob = p_c1_corner * p_win_given_c1_corner + p_c1_edge * p_win_given_c1_edge
    
    print("--- Final Probability Calculation ---")
    print("The maximum chance of winning is the sum of probabilities of these two cases:")
    print(f"P(Win) = P(C1=corner) * P(Win|C1=corner) + P(C1=edge) * P(Win|C1=edge)")
    print(f"P(Win) = {p_c1_corner} * {p_win_given_c1_corner} + {p_c1_edge} * {p_win_given_c1_edge}")
    print(f"P(Win) = {p_c1_corner * p_win_given_c1_corner} + {p_c1_edge * p_win_given_c1_edge}")
    print(f"P(Win) = {total_win_prob}")

    print(f"\nThe maximum chance of winning you can give yourself is {total_win_prob}.")

solve_tic_tac_toe_probability()
<<<8/9>>>