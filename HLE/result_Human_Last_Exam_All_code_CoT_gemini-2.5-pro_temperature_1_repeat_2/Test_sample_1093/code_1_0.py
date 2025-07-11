import random

def run_one_gauss_game():
    """
    Simulates one round of the game as described by Gauss.

    In each game, we sum random numbers between 0 and 1 until the sum
    is larger than 1. We then record the count of numbers used.

    Returns:
        int: The number of random numbers summed to exceed 1.
    """
    current_sum = 0.0
    count = 0
    # Steps 1 & 2: Repeatedly add random numbers until the sum > 1
    while current_sum <= 1:
        random_number = random.random()
        current_sum += random_number
        count += 1
    
    # Step 3: Record the number (the count)
    return count

def main():
    """
    Main function to run the simulation multiple times and report the average.
    """
    # A large number of simulations gives a more accurate average.
    num_simulations = 2000000 
    
    total_recorded_numbers_sum = 0
    
    # Step 4: Reset the sum and repeat the game
    print(f"Running {num_simulations} simulations of the game...")
    for _ in range(num_simulations):
        total_recorded_numbers_sum += run_one_gauss_game()
    
    # Calculate the average of the recorded numbers
    average_value = total_recorded_numbers_sum / num_simulations
    
    # Print the final results and the equation used
    print("\n--- Simulation Results ---")
    print(f"The final equation for the average is:")
    print(f"Average = Total Sum of Recorded Numbers / Number of Simulations")
    print(f"Average = {total_recorded_numbers_sum} / {num_simulations}")
    print(f"The calculated average is: {average_value}")
    print("\nThis value converges to Euler's number, e (approximately 2.71828).")

if __name__ == "__main__":
    main()
