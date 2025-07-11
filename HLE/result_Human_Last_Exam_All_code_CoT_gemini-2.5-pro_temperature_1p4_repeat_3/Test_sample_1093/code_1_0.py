import random

def simulate_one_game():
    """
    Simulates one round of Gauss's game.

    This function repeatedly generates a random number between 0 and 1 and adds it
    to a running sum. It counts how many numbers are needed for the sum to
    become larger than 1.

    Returns:
        int: The number of random numbers (throws) required.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1.
        measurement = random.random()
        
        # Add the measurement to the sum and count the throw.
        current_sum += measurement
        num_throws += 1
    
    # Step 3: Return the recorded number of throws for this game.
    return num_throws

def main():
    """
    Runs the simulation many times and calculates the average number of throws.
    """
    # Step 4: We run a large number of simulations to get a good approximation.
    num_simulations = 2000000
    total_throws_recorded = 0

    print("Running simulation... (This may take a moment)")
    
    for _ in range(num_simulations):
        result = simulate_one_game()
        total_throws_recorded += result
        
    # Calculate the final average.
    average_result = total_throws_recorded / num_simulations

    print("Simulation complete.")
    print("\nThe average is calculated as: Total Throws / Number of Games")
    print(f"{total_throws_recorded} / {num_simulations} = {average_result}")

if __name__ == "__main__":
    main()
