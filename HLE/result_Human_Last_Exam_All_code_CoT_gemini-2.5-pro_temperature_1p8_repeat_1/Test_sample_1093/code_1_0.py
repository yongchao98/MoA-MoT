import random
import math

def run_gauss_game_simulation():
    """
    Simulates Gauss's game for a large number of trials to find the average
    number of throws needed for the sum to exceed 1.
    """
    # Set the number of trials for the simulation. More trials lead to a more accurate result.
    num_trials = 2000000

    # Set a seed for the random number generator for reproducible results.
    random.seed(42)

    # This variable will store the sum of all recorded numbers from each trial.
    total_recorded_numbers = 0

    # Main loop to run the simulation num_trials times.
    for _ in range(num_trials):
        # Step 4: Reset the sum for a new game.
        current_sum = 0.0
        recorded_number = 0
        
        # Step 1-3: Add random numbers until the sum is larger than 1.
        while current_sum <= 1:
            # Generate a random number between 0 and 1.
            random_val = random.random()
            current_sum += random_val
            recorded_number += 1
        
        # Add the recorded number (the count of throws) to the total sum.
        total_recorded_numbers += recorded_number

    # Calculate the average of the recorded numbers over all trials.
    average = total_recorded_numbers / num_trials

    # Print the results of the simulation.
    # The final equation shows the numbers used to calculate the average.
    print("Simulating Gauss's game...")
    print(f"Total sum of recorded numbers: {total_recorded_numbers}")
    print(f"Total number of trials: {num_trials}")
    print(f"Final Equation: {total_recorded_numbers} / {num_trials} = {average}")
    print(f"\nThe simulation shows the average converges to: {average}")
    print(f"The theoretical value is e = {math.e}")
    print(f"Difference from e: {abs(average - math.e)}")

# Run the simulation.
run_gauss_game_simulation()