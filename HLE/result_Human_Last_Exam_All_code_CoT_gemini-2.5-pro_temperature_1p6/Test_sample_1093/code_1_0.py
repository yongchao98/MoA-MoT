import random

def run_simulation():
    """
    This function simulates Gauss's game to find the average number of random
    measurements (numbers between 0 and 1) that must be summed to exceed 1.

    The theoretical average is Euler's number, e (approx. 2.71828).
    This simulation will demonstrate this by running the game many times and
    calculating the average.
    """
    # The number of times we repeat the entire game to get a stable average.
    num_simulations = 1000000
    
    # This list will store the recorded number (n) from each game.
    recorded_numbers = []

    print(f"Running {num_simulations} simulations of the game...")

    # Main loop to run the simulations.
    for _ in range(num_simulations):
        current_sum = 0.0
        num_throws = 0
        # Loop for a single game.
        while current_sum <= 1:
            # Step 1 & 2: Generate a random number between 0 and 1.
            measurement = random.random()
            current_sum += measurement
            num_throws += 1
        # Step 3: Record the number.
        recorded_numbers.append(num_throws)

    # Calculate the final average.
    total_sum_of_throws = sum(recorded_numbers)
    number_of_games = len(recorded_numbers)
    average_result = total_sum_of_throws / number_of_games

    # Print the final equation and the result.
    print("\nCalculation of the average:")
    print(f"The sum of all recorded numbers is: {total_sum_of_throws}")
    print(f"The number of games played is: {number_of_games}")
    print(f"Average = {total_sum_of_throws} / {number_of_games} = {average_result}")
    print("\nThis average is an approximation of Euler's number, e.")

run_simulation()