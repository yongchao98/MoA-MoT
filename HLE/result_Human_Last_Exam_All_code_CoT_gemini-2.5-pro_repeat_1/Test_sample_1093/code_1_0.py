import random

def run_one_gauss_game():
    """
    Simulates one round of Gauss's game.
    In each step, a random number between 0 and 1 is generated and added to a sum.
    The function returns the number of steps it took for the sum to exceed 1.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1
        measurement = random.random()
        # Step 3: Add to the sum
        current_sum += measurement
        num_throws += 1
    # Record the number of throws
    return num_throws

def main():
    """
    Runs the simulation multiple times and calculates the average number of throws.
    """
    # Step 4: Repeat the game many times
    num_simulations = 200000
    total_throws_recorded = 0
    
    for _ in range(num_simulations):
        total_throws_recorded += run_one_gauss_game()

    # Calculate the average
    average_throws = total_throws_recorded / num_simulations
    
    print(f"Simulation ran for {num_simulations} games.")
    print("The final average is calculated as:")
    # Output the final equation with the numbers
    print(f"{total_throws_recorded} / {num_simulations} = {average_throws}")

if __name__ == "__main__":
    main()
