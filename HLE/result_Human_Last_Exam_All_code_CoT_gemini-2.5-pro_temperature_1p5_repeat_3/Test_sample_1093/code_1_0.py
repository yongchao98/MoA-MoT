import random

def run_gauss_game_simulation(num_simulations):
    """
    Simulates J.C.F. Gauss's game a specified number of times to find the
    average number of throws needed for the sum to exceed 1.

    Args:
        num_simulations (int): The total number of games to play.
    """
    # A seed ensures the "random" numbers are the same each time the code is run,
    # making the result reproducible.
    random.seed(42)

    total_throws_across_all_games = 0

    # Step 4: Repeat the game for the specified number of simulations
    for _ in range(num_simulations):

        # --- This block simulates a single game ---
        current_sum = 0.0
        throws_in_this_game = 0
        while current_sum <= 1:
            # Step 1-2: Generate a random number and add it to the sum
            random_value = random.random()
            current_sum += random_value
            throws_in_this_game += 1
        # --- End of single game simulation ---

        # Step 3: Record the number by adding it to our total
        total_throws_across_all_games += throws_in_this_game

    # Calculate the final average
    average_result = total_throws_across_all_games / num_simulations

    # Print the final equation with all its components
    print(f"Total throws from all games: {total_throws_across_all_games}")
    print(f"Total number of games played: {num_simulations}")
    print(f"Final Equation: {total_throws_across_all_games} / {num_simulations} = {average_result}")
    print("\nThis average converges to Euler's number (e).")


if __name__ == '__main__':
    # We will run 1,000,000 simulations for a good approximation.
    num_games_to_play = 1000000
    run_gauss_game_simulation(num_games_to_play)
