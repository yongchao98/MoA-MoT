import random
import math

def run_gauss_game_trial():
    """
    Simulates one round of Gauss's game.

    In this game, we repeatedly draw random numbers between 0 and 1 and add them
    together until the sum exceeds 1.

    Returns:
        int: The number of draws required for the sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        measurement = random.random()
        
        # Add the measurement to the sum
        current_sum += measurement
        num_throws += 1
    
    # Return the recorded number for this trial
    return num_throws

def main():
    """
    Runs a simulation of Gauss's game to find the average number of throws.
    """
    # A large number of trials will give a result closer to the theoretical average
    num_trials = 1_000_000
    
    print(f"Running a simulation of Gauss's game with {num_trials:,} trials...")
    
    total_throws = 0
    # Step 4: Reset and repeat the game for the specified number of trials
    for _ in range(num_trials):
        total_throws += run_gauss_game_trial()
        
    # Calculate the average of the recorded numbers
    simulated_average = total_throws / num_trials
    
    # The theoretical value is Euler's number, e
    theoretical_average = math.e

    print("\n--- Simulation Results ---")
    print(f"The average of the recorded numbers from the simulation is: {simulated_average}")
    print(f"The theoretical average is Euler's number, e:             {theoretical_average}")
    
    # Show the difference to see how close the simulation is
    difference = abs(simulated_average - theoretical_average)
    print(f"The absolute difference is:                               {difference}")

if __name__ == "__main__":
    main()
