def solve_tichu_puzzle():
    """
    Calculates the maximal possible value of X-Y in a Tichu round
    without a double victory for the winning team.
    """
    print("Let's determine the maximal score difference (X - Y) in a Tichu round.")
    print("X is the winning team's score, and Y is the losing team's score.")
    print("-" * 30)
    
    # Step 1: Maximize the winning team's (Team A) score.
    # To maximize card points, Team A takes all point cards.
    card_points_A = 100
    # To maximize call points, one player from Team A successfully calls Grand Tichu.
    call_points_A = 200
    
    # Calculate Team A's total score (X)
    X = card_points_A + call_points_A
    
    print("To maximize the difference, we maximize the winning team's score (X):")
    print(f"  - Winning team's card points: {card_points_A} (capturing all Kings, Tens, Fives, and the Dragon)")
    print(f"  - Winning team's call points: {call_points_A} (from one successful Grand Tichu call)")
    print(f"Total Winning Score (X) = {card_points_A} + {call_points_A} = {X}")
    print("-" * 30)

    # Step 2: Minimize the losing team's (Team B) score.
    # Since Team A took all card points, Team B has none.
    card_points_B = 0
    # To minimize call points, both players on Team B fail a Grand Tichu call.
    call_points_B = -200 + -200
    
    # Calculate Team B's total score (Y)
    Y = card_points_B + call_points_B
    
    print("Next, we minimize the losing team's score (Y):")
    print(f"  - Losing team's card points: {card_points_B}")
    print(f"  - Losing team's call points: {call_points_B} (from two failed Grand Tichu calls)")
    print(f"Total Losing Score (Y) = {card_points_B} + ({call_points_B}) = {Y}")
    print("-" * 30)

    # Step 3: Calculate the final difference, X - Y.
    max_diff = X - Y
    
    print("Finally, we calculate the maximal value of X - Y.")
    print(f"The equation is: X - Y = ({card_points_A} + {call_points_A}) - ({card_points_B} + {call_points_B})")
    print(f"Final Calculation: X - Y = {X} - ({Y}) = {max_diff}")
    
solve_tichu_puzzle()
<<<700>>>