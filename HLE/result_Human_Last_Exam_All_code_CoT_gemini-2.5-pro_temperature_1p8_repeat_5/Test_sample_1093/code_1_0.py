import random
import math

def run_simulation():
    """
    This function simulates Gauss's game to find the average number of "throws"
    needed for the sum of random numbers (from 0 to 1) to exceed 1.

    The theoretical average is Euler's number, e. This simulation demonstrates
    this principle by running the game many times and calculating the average.
    """
    # The number of times we will repeat the game to get a good average.
    num_simulations = 2000000
    total_throws_recorded = 0

    print(f"Running the simulation {num_simulations:,} times...")

    # Loop for each game simulation.
    for _ in range(num_simulations):
        
        # Reset the sum and throw count for the new game.
        current_sum = 0.0
        throw_count = 0
        
        # Keep throwing until the sum is greater than 1.
        while current_sum <= 1:
            # A "throw" is a random number between 0 and 1.
            throw = random.random()
            current_sum += throw
            throw_count += 1
            
        # Add the number of throws from this game to our total.
        total_throws_recorded += throw_count

    # Calculate the final average.
    average_throws = total_throws_recorded / num_simulations

    # Print the equation used for the final calculation.
    print("\nCalculation of the average:")
    print(f"Total Throws / Number of Games = Average Throws")
    print(f"{total_throws_recorded} / {num_simulations} = {average_throws}")

    # For comparison, print the theoretical value (Euler's number).
    print(f"\nThe simulation result is very close to the theoretical value, e ≈ {math.e}")

if __name__ == "__main__":
    run_simulation()