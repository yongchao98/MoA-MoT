import random

def run_gauss_game_simulation(num_trials):
    """
    Simulates the game described by J.C.F. Gauss to find the average
    number of steps required for the sum of random numbers to exceed 1.

    Args:
        num_trials (int): The number of times to repeat the game for the simulation.
    """
    total_throws_recorded = 0

    # Loop through the number of trials to simulate the game multiple times
    for _ in range(num_trials):
        sum_of_measurements = 0.0
        throws_in_this_game = 0
        while sum_of_measurements <= 1:
            # Step 1 & 2: Generate a random number between 0 and 1
            sum_of_measurements += random.random()
            throws_in_this_game += 1
        # Step 3: Record the number of throws
        total_throws_recorded += throws_in_this_game

    # Calculate the average
    average_throws = total_throws_recorded / num_trials

    print(f"The game was simulated {num_trials:,} times.")
    print("The final equation to calculate the average is:")
    # Output each number in the final equation
    print(f"Average = Total Recorded Throws / Number of Games")
    print(f"Average = {total_throws_recorded} / {num_trials} = {average_throws}")
    print("\nAs you can see, the average converges to Euler's number, e (approximately 2.71828).")


if __name__ == '__main__':
    # Step 4: Run the simulation a large number of times for an accurate result
    number_of_simulations = 2000000
    run_gauss_game_simulation(number_of_simulations)
