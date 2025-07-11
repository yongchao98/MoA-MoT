from fractions import Fraction

def calculate_tic_tac_toe_win_prob():
    """
    Calculates and prints the maximum win probability in Tic Tac Toe
    against a random opponent, assuming a tie is a loss.
    """

    print("Analyzing the optimal strategy for Tic Tac Toe against a random opponent.")
    print("A tie is considered a loss. We must maximize our win probability.")
    print("-" * 60)

    # --- Strategy 1: Start in a Corner ---
    print("Strategy 1: Player starts in a corner.")

    # Probability of O playing center on its first move
    prob_o_plays_center = Fraction(1, 8)
    # Probability of winning if O plays center (requires deeper analysis not shown here)
    win_prob_if_o_center = Fraction(5, 6)

    # Probability of O NOT playing center on its first move
    prob_o_not_center = Fraction(7, 8)
    # Probability of winning if O does not play center (player can force a win)
    win_prob_if_o_not_center = Fraction(1, 1)

    # Total win probability for starting in a corner
    total_win_prob_corner = (prob_o_plays_center * win_prob_if_o_center) + \
                            (prob_o_not_center * win_prob_if_o_not_center)

    print(f"The computer has 8 opening moves.")
    print(f"  - P(O plays center) = {prob_o_plays_center}")
    print(f"    - In this case, our optimal play leads to a win probability of {win_prob_if_o_center}.")
    print(f"  - P(O does not play center) = {prob_o_not_center}")
    print(f"    - In this case, we can force a win, so the win probability is {win_prob_if_o_not_center}.")
    
    print("\nTotal win probability for starting in a corner:")
    # Print the equation with all numbers
    print(f"P(win) = ({prob_o_plays_center}) * ({win_prob_if_o_center}) + ({prob_o_not_center}) * ({win_prob_if_o_not_center})")
    print(f"P(win) = {prob_o_plays_center * win_prob_if_o_center} + {prob_o_not_center * win_prob_if_o_not_center}")
    print(f"P(win) = {total_win_prob_corner}")
    print("-" * 60)

    # --- Strategy 2: Start in the Center ---
    print("Strategy 2: Player starts in the center.")

    # Probability of O playing a corner on its first move
    prob_o_plays_corner = Fraction(4, 8)
    # Win probability if O plays a corner
    win_prob_if_o_corner = Fraction(3, 4)

    # Probability of O playing an edge on its first move
    prob_o_plays_edge = Fraction(4, 8)
    # Win probability if O plays an edge (player can force a win)
    win_prob_if_o_edge = Fraction(1, 1)

    # Total win probability for starting in the center
    total_win_prob_center = (prob_o_plays_corner * win_prob_if_o_corner) + \
                            (prob_o_plays_edge * win_prob_if_o_edge)

    print(f"The computer has 8 opening moves.")
    print(f"  - P(O plays a corner) = {prob_o_plays_corner}")
    print(f"    - In this case, our optimal play leads to a win probability of {win_prob_if_o_corner}.")
    print(f"  - P(O plays an edge) = {prob_o_plays_edge}")
    print(f"    - In this case, we can force a win, so the win probability is {win_prob_if_o_edge}.")
    
    print("\nTotal win probability for starting in the center:")
    # Print the equation with all numbers
    print(f"P(win) = ({prob_o_plays_corner}) * ({win_prob_if_o_corner}) + ({prob_o_plays_edge}) * ({win_prob_if_o_edge})")
    print(f"P(win) = {prob_o_plays_corner * win_prob_if_o_corner} + {prob_o_plays_edge * win_prob_if_o_edge}")
    print(f"P(win) = {total_win_prob_center}")
    print("-" * 60)

    # --- Comparison and Conclusion ---
    print("Comparing the two strategies:")
    print(f"  - Starting Corner P(win) = {total_win_prob_corner}")
    print(f"  - Starting Center P(win) = {total_win_prob_center}")

    if total_win_prob_corner > total_win_prob_center:
        max_prob = total_win_prob_corner
        best_strategy = "starting in a corner"
    else:
        max_prob = total_win_prob_center
        best_strategy = "starting in the center"

    print(f"\nThe maximum chance of winning is achieved by {best_strategy}, which gives a probability of {max_prob}.")


calculate_tic_tac_toe_win_prob()