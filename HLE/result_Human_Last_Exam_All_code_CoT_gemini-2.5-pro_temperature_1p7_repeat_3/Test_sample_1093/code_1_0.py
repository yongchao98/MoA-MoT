import random

def run_simulation():
    """
    This function simulates the game described by J.C.F. Gauss.

    The game involves:
    1. Generating a random number between 0 and 1 (simulating the needle throw).
    2. Adding these numbers until the sum is greater than 1.
    3. Recording the count of numbers needed.

    The function simulates this process many times and calculates the average
    of the recorded counts. This average is known to converge to Euler's number, e.
    """
    # The number of games to simulate. A larger number gives a more accurate average.
    num_games = 1000000

    total_throws_recorded = 0

    # Loop to play the game num_games times
    for _ in range(num_games):
        current_sum = 0.0
        throws_for_this_game = 0
        while current_sum <= 1:
            # Generate a random number between 0 and 1
            random_value = random.random()
            current_sum += random_value
            throws_for_this_game += 1
        
        # Add the number of throws from this game to the total
        total_throws_recorded += throws_for_this_game

    # Calculate the average number of throws
    average_throws = total_throws_recorded / num_games

    print(f"The simulation was run for {num_games} games.")
    print("The theoretical average is Euler's number, e ≈ 2.71828...")
    print("\n--- Simulation Results ---")
    # Here we output the numbers used in the final calculation
    print(f"Total throws recorded: {total_throws_recorded}")
    print(f"Number of games: {num_games}")
    print(f"Average number of throws = {total_throws_recorded} / {num_games} = {average_throws}")

# Run the simulation
run_simulation()