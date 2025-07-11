import random
import math

def run_gauss_game():
    """
    Simulates a single game, summing random numbers until the total exceeds 1.
    Returns the count of numbers (throws) required.
    """
    sum_of_throws = 0.0
    number_of_throws = 0
    while sum_of_throws <= 1:
        # A random number between 0 and 1, as per the problem description.
        throw = random.random()
        sum_of_throws += throw
        number_of_throws += 1
    return number_of_throws

# Use a large number of games to get a good approximation of the true average.
num_simulations = 2_000_000
total_throws_recorded = 0

# Run the simulation many times.
for _ in range(num_simulations):
    total_throws_recorded += run_gauss_game()

# Calculate the average number of throws.
average_throws = total_throws_recorded / num_simulations

# The theoretical answer is Euler's number, e.
# We will print our simulated result to show how it converges to e.
print("This simulation estimates the average number of throws needed for their sum to exceed 1.")
print(f"Number of simulated games: {num_simulations}")
print(f"Total throws across all games: {total_throws_recorded}")
print("\nThe average is calculated as: Total Throws / Number of Games")
print(f"Final Equation: {total_throws_recorded} / {num_simulations} = {average_throws}")

print(f"\nThe simulation result of {average_throws} is an approximation of the theoretical value.")
print(f"The exact value the average converges to is Euler's number, e ≈ {math.e}")
