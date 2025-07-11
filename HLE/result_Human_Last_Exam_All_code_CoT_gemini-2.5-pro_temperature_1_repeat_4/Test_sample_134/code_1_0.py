import math

def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # 1. Game Parameters
    n1 = 136  # Number of 1-euro coins
    n2 = 87   # Number of 2-euro coins
    N = n1 + n2  # Total number of coins

    # Total turns is N, which is odd.
    # Player 1 takes (N+1)/2 turns, Player 2 takes (N-1)/2 turns.
    p1_turns = (N + 1) // 2
    p2_turns = (N - 1) // 2

    print("Step 1: Understand the game structure.")
    print(f"There are a total of {N} coins ({n1} x 1-euro, {n2} x 2-euro).")
    print(f"The total number of turns is {N}, which is an odd number.")
    print(f"This means Player 1 gets to pick {p1_turns} coins.")
    print(f"Player 2 gets to pick {p2_turns} coins.")
    print("Player 1 always picks one more coin than Player 2.")
    print("-" * 40)

    # 2. The Key Insight: Coloring the Positions
    n_odd = p1_turns
    n_even = p2_turns
    
    print("Step 2: The 'Position Coloring' Strategy.")
    print("Let's label coin positions as 'odd' (1, 3, 5, ...) and 'even' (2, 4, 6, ...).")
    print(f"There are {n_odd} odd positions and {n_even} even positions in the line of coins.")
    print("A key observation is that on Player 1's turn, both available coins at the ends are on positions of the SAME parity (e.g., pos 1 and pos 223 are both odd).")
    print("On Player 2's turn, the available coins are on positions of DIFFERENT parity (one odd, one even).")
    print("-" * 40)

    # 3. Player 2's Control
    print("Step 3: Player 2's control over the game.")
    print("Because Player 2 always chooses between an odd-position coin and an even-position coin, they have significant control.")
    print("Let S_odd be the sum of values of all coins in odd positions, and S_even for even positions.")
    print("Player 2 can adopt an optimal strategy based on these sums for any given arrangement of coins:")
    print(" - If S_even > S_odd, P2 can guarantee a win by always choosing the even-position coin. This forces P1 to take from the less valuable odd-position pool.")
    print(" - If S_odd >= S_even, P2's best strategy is to always choose the odd-position coin. This limits P1 to taking from the even-position pool (except for P1's very first pick).")
    print("-" * 40)

    # 4. Analyzing the score difference
    print("Step 4: Analyzing the score difference based on the coin arrangement.")
    print("The arrangement is random. Let 'x' be the number of 2-euro coins that fall in the {n_odd} odd positions.")
    print("We can express the sums S_odd and S_even in terms of x.")
    # S_odd = 1*(num_1s_in_odd) + 2*(num_2s_in_odd) = 1*(n_odd - x) + 2*x = n_odd + x
    # Total Value = 1*136 + 2*87 = 136 + 174 = 310
    # S_even = Total Value - S_odd = 310 - (n_odd + x)
    # S_odd - S_even = (n_odd + x) - (310 - (n_odd + x)) = 2*n_odd + 2*x - 310
    # S_odd - S_even = 2*{n_odd} + 2*x - {n1 + 2*n2} = 2*112 + 2x - 310 = 224 + 2x - 310 = 2x - 86
    print(f"The difference in total value is: S_odd - S_even = 2*x - 86.")
    print("-" * 40)
    
    # 5. Winning Conditions based on x
    print("Step 5: Winning conditions based on x.")
    print(f"Case A: S_even > S_odd  ==>  2*x - 86 < 0  ==>  x < 43.")
    print("   In this case, Player 2 employs their strategy and wins.")
    print(f"Case B: S_odd = S_even  ==>  2*x - 86 = 0  ==>  x = 43.")
    print("   Here, Player 1's score will be S_even + c_max and Player 2's will be S_odd - c_max (where c_max is the value of the coin P1 picks first).")
    print("   Since S_odd = S_even, Player 1 wins by a score difference of 2 * c_max.")
    print(f"Case C: S_odd > S_even  ==>  2*x - 86 > 0  ==>  x > 43.")
    print("   Here again, Score_P1 = S_even + c_max and Score_P2 = S_odd - c_max.")
    print("   Player 2 wins only if Score_P2 > Score_P1, which means (S_odd - S_even) > 2*c_max.")
    print("   Substituting our formula: (2*x - 86) > 2*c_max  ==>  x - 43 > c_max.")
    print("-" * 40)

    # 6. Most Likely Scenario
    # The number of 2-euro coins in odd positions, x, follows a hypergeometric distribution.
    # We find its mode to see the most likely value of x.
    N_pop = N         # Population size = 223 coins
    K_pop = n2        # Number of "successes" in pop = 87 2-euro coins
    n_sample = n_odd  # Sample size = 112 odd positions
    mode = math.floor((n_sample + 1) * (K_pop + 1) / (N_pop + 2))
    
    print("Step 6: Find the most likely arrangement.")
    print("The variable 'x' follows a hypergeometric distribution. We find its most probable value (the mode).")
    print(f"Mode = floor( (n_sample + 1) * (K_pop + 1) / (N_pop + 2) )")
    print(f"Mode = floor( ({n_sample} + 1) * ({K_pop} + 1) / ({N_pop} + 2) ) = floor( {113 * 88} / {225} ) = floor({9944 / 225:.2f}) = {mode}")
    print(f"This means the most likely number of 2-euro coins in odd positions is x = {mode}.")
    print("-" * 40)
    
    # 7. Analyze the outcome for the mode and conclude
    print("Step 7: Analyze the outcome for the most likely scenarios and conclude.")
    print(f"The most probable arrangement of coins is the one where x = {mode}.")
    print(f"For x = 44, we are in Case C (since 44 > 43).")
    print("   Player 2 wins if x - 43 > c_max, which for x=44 becomes 44 - 43 > c_max ==> 1 > c_max.")
    print("   This is impossible, as the coin value c_max is at least 1. So Player 2 can't win.")
    print("   If c_max = 1, the game is a tie. If c_max = 2, Player 1 wins.")
    print("   Therefore, in the most likely scenario, Player 1 either wins or ties, but never loses.")
    print("\nThe next most likely scenario is x = 43.")
    print("   In this case, as determined in Case B, Player 1 wins for sure.")
    print("\nCONCLUSION:")
    print("Since the most probable random arrangements of coins lead to a situation where Player 1 either wins or ties, it is strategically better to choose to be the first player.")

# Run the analysis
solve_coin_game()