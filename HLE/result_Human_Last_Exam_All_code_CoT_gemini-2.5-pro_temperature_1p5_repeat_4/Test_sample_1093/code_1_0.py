import random

def simulate_gauss_game():
    """
    Simulates the game described by J.C.F. Gauss to find the average
    number of throws needed for the sum of random numbers to exceed 1.
    """
    # Set the number of games (simulations) to run for a good approximation.
    num_simulations = 1_000_000
    
    # This list will store the result (number of throws) of each game.
    recorded_numbers = []

    # Loop for each simulation (game).
    for _ in range(num_simulations):
        # Step 4: Reset the sum and repeat for a new game.
        current_sum = 0.0
        num_throws = 0
        
        # Step 3: Add measurements until the sum is larger than 1.
        while current_sum <= 1:
            # Step 1-2: Generate a random number between 0 and 1.
            # This corresponds to the normalized angle.
            random_value = random.random()
            
            current_sum += random_value
            num_throws += 1
        
        # Record the number of throws for this completed game.
        recorded_numbers.append(num_throws)

    # Calculate the total sum of throws across all games.
    total_throws_sum = sum(recorded_numbers)
    
    # Calculate the average.
    average_result = total_throws_sum / num_simulations

    # Print the final calculation and the result, as requested.
    print("The average is calculated as the sum of all recorded throw counts divided by the total number of games.")
    print(f"Average = {total_throws_sum} / {num_simulations}")
    print(f"The average of the recorded numbers converges to: {average_result}")

if __name__ == '__main__':
    simulate_gauss_game()