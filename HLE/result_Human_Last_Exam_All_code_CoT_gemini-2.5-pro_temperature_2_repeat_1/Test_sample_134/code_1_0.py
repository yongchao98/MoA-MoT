import collections

def solve_coin_game():
    """
    Solves the coin game problem by simulating an example and explaining the logic.
    """

    # Step 1: Explain the core logic based on the number of turns.
    total_coins = 136 + 87
    player1_coins = (total_coins + 1) // 2
    player2_coins = (total_coins - 1) // 2
    total_value = 136 * 1 + 87 * 2

    # Explanation text
    print("Thinking Process:")
    print(f"Total number of coins is 136 + 87 = {total_coins}, which is an odd number.")
    print(f"This means Player 1 gets to pick {player1_coins} coins, while Player 2 gets to pick {player2_coins} coins.")
    print("This inherent advantage of getting an extra coin means Player 1 has a guaranteed non-losing strategy.")
    print("To demonstrate this, let's simulate the game with one specific arrangement:")
    print("The arrangement is 44 '2-euro' coins, then 136 '1-euro' coins, then 43 '2-euro' coins.\n")

    # Step 2: Simulate the game with the chosen arrangement.
    # The deque is efficient for removing items from both ends.
    coins = collections.deque([2] * 44 + [1] * 136 + [2] * 43)
    
    player1_score = 0
    player2_score = 0
    
    for turn in range(total_coins):
        is_player1_turn = (turn % 2 == 0)

        # The optimal move is to take the most valuable coin available.
        left_coin = coins[0]
        right_coin = coins[-1]
        
        # Decide which coin to pick
        if left_coin > right_coin:
            picked_coin = coins.popleft()
        elif right_coin > left_coin:
            picked_coin = coins.pop()
        else:
            # If both ends have the same value, it doesn't matter which is picked.
            picked_coin = coins.popleft()
            
        if is_player1_turn:
            player1_score += picked_coin
        else:
            player2_score += picked_coin

    # Step 3: Print the results of the simulation.
    print("Simulation Result:")
    print(f"Player 1's final score: {player1_score}")
    print(f"Player 2's final score: {player2_score}")
    print("\nFinal Equation:")
    print(f"{player1_score} + {player2_score} = {player1_score + player2_score}")

solve_coin_game()