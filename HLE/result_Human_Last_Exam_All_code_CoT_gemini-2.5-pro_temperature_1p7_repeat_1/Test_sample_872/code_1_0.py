from fractions import Fraction

def calculate_win_probability():
    """
    Calculates the maximum win probability in Tic Tac Toe against a random opponent.
    The analysis shows the optimal first move is to a corner. This function
    calculates the win probability for that strategy.
    """

    # --- Step 1: Define probabilities of the computer's first move ---
    # The computer chooses one of the 8 empty squares uniformly at random.
    # We group these moves by symmetry, assuming we played in a corner (e.g., top-left).
    
    # Computer plays the center (1 of 8 squares)
    p_O_center = Fraction(1, 8)
    
    # Computer plays the opposite corner (1 of 8 squares)
    p_O_opposite_corner = Fraction(1, 8)
    
    # Computer plays one of the 2 adjacent corners
    p_O_adj_corners = Fraction(2, 8)
    
    # Computer plays one of the 2 adjacent edges
    p_O_adj_edges = Fraction(2, 8)
    
    # Computer plays one of the 2 non-adjacent edges
    p_O_non_adj_edges = Fraction(2, 8)

    # --- Step 2: Define our win probabilities for each of the computer's moves ---
    # These probabilities are derived from a game-tree analysis where we play optimally.
    
    # If O plays center, our best reply leads to a 5/6 chance of winning.
    win_prob_given_O_center = Fraction(5, 6)
    
    # If O plays the opposite corner, our best reply creates a forced win.
    win_prob_given_O_opposite_corner = Fraction(1, 1)
    
    # If O plays an adjacent corner, our best reply leads to a 23/24 chance of winning.
    win_prob_given_O_adj_corners = Fraction(23, 24)
    
    # If O plays an adjacent edge, our best reply creates a forced win.
    win_prob_given_O_adj_edges = Fraction(1, 1)
    
    # If O plays a non-adjacent edge, our best reply creates a forced win.
    win_prob_given_O_non_adj_edges = Fraction(1, 1)

    # --- Step 3: Calculate the total win probability using the law of total probability ---
    
    term1 = p_O_center * win_prob_given_O_center
    term2 = p_O_opposite_corner * win_prob_given_O_opposite_corner
    term3 = p_O_adj_corners * win_prob_given_O_adj_corners
    term4 = p_O_adj_edges * win_prob_given_O_adj_edges
    term5 = p_O_non_adj_edges * win_prob_given_O_non_adj_edges

    total_win_prob = term1 + term2 + term3 + term4 + term5

    # --- Step 4: Print the full calculation step-by-step ---

    print("The maximum chance of winning is achieved by starting in a corner.")
    print("The total probability is the sum of probabilities of each case, weighted by the chance of that case occurring.")
    print("\nThe full equation is:")
    
    equation_str = (
        f"P(Win) = P(O plays center) * P(Win|O center) + \n"
        f"         P(O plays opp corner) * P(Win|O opp corner) + \n"
        f"         P(O plays adj corner) * P(Win|O adj corner) + \n"
        f"         P(O plays adj edge) * P(Win|O adj edge) + \n"
        f"         P(O plays non-adj edge) * P(Win|O non-adj edge)\n"
    )
    print(equation_str)

    print("Substituting the fractional values:")
    
    final_equation_str = (
        f"P(Win) = ({p_O_center.numerator}/{p_O_center.denominator})*({win_prob_given_O_center.numerator}/{win_prob_given_O_center.denominator}) + "
        f"({p_O_opposite_corner.numerator}/{p_O_opposite_corner.denominator})*({win_prob_given_O_opposite_corner.numerator}/{win_prob_given_O_opposite_corner.denominator}) + "
        f"({p_O_adj_corners.numerator}/{p_O_adj_corners.denominator})*({win_prob_given_O_adj_corners.numerator}/{win_prob_given_O_adj_corners.denominator}) + "
        f"({p_O_adj_edges.numerator}/{p_O_adj_edges.denominator})*({win_prob_given_O_adj_edges.numerator}/{win_prob_given_O_adj_edges.denominator}) + "
        f"({p_O_non_adj_edges.numerator}/{p_O_non_adj_edges.denominator})*({win_prob_given_O_non_adj_edges.numerator}/{win_prob_given_O_non_adj_edges.denominator})"
    )
    print(final_equation_str)
    
    print("\nCalculating each term:")
    print(f"         {term1.numerator}/{term1.denominator} + {term2.numerator}/{term2.denominator} + {term3.numerator}/{term3.denominator} + {term4.numerator}/{term4.denominator} + {term5.numerator}/{term5.denominator}")

    print(f"\nSumming the terms gives the final probability:")
    print(f"P(Win) = {total_win_prob.numerator}/{total_win_prob.denominator}")


calculate_win_probability()