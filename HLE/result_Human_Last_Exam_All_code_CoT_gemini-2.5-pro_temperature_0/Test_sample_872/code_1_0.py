from fractions import Fraction

def solve_tic_tac_toe_probability():
    """
    Calculates the maximum probability of winning Tic Tac Toe against a random computer.
    """

    # The optimal strategy is to start in the center.
    # The computer has 8 squares left. 4 are corners, 4 are edges.
    prob_comp_plays_corner = Fraction(4, 8)
    prob_comp_plays_edge = Fraction(4, 8)

    print("My optimal first move is to take the center square.")
    print(f"The computer can then play a corner (probability {prob_comp_plays_corner}) or an edge (probability {prob_comp_plays_edge}).\n")

    # Case 1: Computer plays an edge.
    # My response creates a threat. The computer has 6 squares to choose from to block.
    prob_comp_fails_first_block_edge = Fraction(5, 6)
    # If the computer blocks, my next move creates a fork, guaranteeing a win.
    prob_win_if_comp_blocks_edge = Fraction(1, 1)
    prob_win_if_comp_plays_edge = prob_comp_fails_first_block_edge * 1 + (1 - prob_comp_fails_first_block_edge) * prob_win_if_comp_blocks_edge
    
    print("--- Analysis if Computer Plays an Edge ---")
    print(f"If the computer plays an edge, I can create a threat.")
    print(f"The computer has a {prob_comp_fails_first_block_edge} chance of failing to block, in which case I win.")
    print(f"If it blocks, my next move creates a fork, guaranteeing a win (probability {prob_win_if_comp_blocks_edge}).")
    print(f"So, P(Win | Computer plays Edge) = {prob_comp_fails_first_block_edge} * 1 + {1-prob_comp_fails_first_block_edge} * {prob_win_if_comp_blocks_edge} = {prob_win_if_comp_plays_edge}\n")

    # Case 2: Computer plays a corner.
    # My response creates a threat. Computer has 6 squares to block.
    prob_comp_fails_o2_block = Fraction(5, 6)
    prob_comp_makes_o2_block = Fraction(1, 6)

    # If O2 blocks, I make a new threat. O3 has 4 squares to block.
    prob_comp_fails_o3_block = Fraction(3, 4)
    prob_comp_makes_o3_block = Fraction(1, 4)

    # If O3 blocks, I make a new threat. O4 has 2 squares to block.
    prob_comp_fails_o4_block = Fraction(1, 2)
    prob_comp_makes_o4_block = Fraction(1, 2)
    
    # If O4 blocks, the game is a tie (win probability = 0).
    prob_win_after_o4_block = Fraction(0, 1)

    # Calculate backwards
    prob_win_after_o3_block = prob_comp_fails_o4_block * 1 + prob_comp_makes_o4_block * prob_win_after_o4_block
    prob_win_after_o2_block = prob_comp_fails_o3_block * 1 + prob_comp_makes_o3_block * prob_win_after_o3_block
    prob_win_if_comp_plays_corner = prob_comp_fails_o2_block * 1 + prob_comp_makes_o2_block * prob_win_after_o2_block

    print("--- Analysis if Computer Plays a Corner ---")
    print("The game tree is more complex. We calculate the win probability by considering each of the computer's random moves.")
    print(f"P(Win after O3 blocks) = {prob_comp_fails_o4_block} * 1 + {prob_comp_makes_o4_block} * 0 = {prob_win_after_o3_block}")
    print(f"P(Win after O2 blocks) = {prob_comp_fails_o3_block} * 1 + {prob_comp_makes_o3_block} * {prob_win_after_o3_block} = {prob_win_after_o2_block}")
    print(f"P(Win | Computer plays Corner) = {prob_comp_fails_o2_block} * 1 + {prob_comp_makes_o2_block} * {prob_win_after_o2_block} = {prob_win_if_comp_plays_corner}\n")

    # Final Calculation
    total_win_prob = prob_win_if_comp_plays_corner * prob_comp_plays_corner + prob_win_if_comp_plays_edge * prob_comp_plays_edge

    print("--- Final Calculation ---")
    print("P(Win) = P(Win | Corner) * P(Corner) + P(Win | Edge) * P(Edge)")
    print(f"P(Win) = {prob_win_if_comp_plays_corner} * {prob_comp_plays_corner} + {prob_win_if_comp_plays_edge} * {prob_comp_plays_edge}")
    print(f"P(Win) = {prob_win_if_comp_plays_corner * prob_comp_plays_corner} + {prob_win_if_comp_plays_edge * prob_comp_plays_edge}")
    print(f"P(Win) = {total_win_prob}")
    
    print(f"\nThe maximum chance of winning is {total_win_prob.numerator}/{total_win_prob.denominator}.")

solve_tic_tac_toe_probability()