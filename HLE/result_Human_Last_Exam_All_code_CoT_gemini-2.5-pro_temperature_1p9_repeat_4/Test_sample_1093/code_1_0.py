import random

def simulate_gauss_game():
    """
    This script simulates the needle game to find the average number of throws
    needed for the sum of normalized angles to exceed 1.
    """
    # Number of times to repeat the game to get a good average.
    num_simulations = 1_000_000
    
    total_throws_recorded = 0

    # Main loop to run the simulation multiple times.
    for _ in range(num_simulations):
        current_sum = 0.0
        num_throws = 0
        # A single game: keep adding numbers until the sum is > 1.
        while current_sum <= 1:
            # Generate a random number between 0 and 1.
            random_value = random.random()
            current_sum += random_value
            num_throws += 1
        
        # Add the recorded number of throws from this game to the total.
        total_throws_recorded += num_throws

    # Calculate the average.
    average_throws = total_throws_recorded / num_simulations

    print("This script simulates the game to find what the average of the recorded numbers converges to.")
    print(f"Number of games simulated: {num_simulations}")
    
    # Per the instructions, printing the numbers used in the final equation.
    print("\nFinal Average Calculation:")
    print(f"{total_throws_recorded} (Total Throws) / {num_simulations} (Number of Games) = {average_throws}")
    
    print("\nAs shown by the simulation, the average converges to Euler's number, e ≈ 2.71828.")

simulate_gauss_game()