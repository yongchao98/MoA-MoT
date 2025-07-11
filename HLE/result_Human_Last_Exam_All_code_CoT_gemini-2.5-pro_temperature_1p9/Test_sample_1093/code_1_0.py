import random
import math

def simulate_one_game():
    """
    Simulates one round of Gauss's game.
    In each round, we sum random numbers (from 0 to 1) until the sum exceeds 1.
    This function returns the number of random numbers it took.
    """
    current_sum = 0.0
    num_throws = 0
    while current_sum <= 1:
        # Step 1 & 2: Generate a random number between 0 and 1 and add to sum
        random_value = random.random()
        current_sum += random_value
        num_throws += 1
    # Step 3: Record the number of throws
    return num_throws

# Step 4: Reset and repeat many times
# A large number of simulations gives a better approximation of the average.
num_simulations = 2000000
total_throws_recorded = 0

# Run the simulation loop
for _ in range(num_simulations):
    total_throws_recorded += simulate_one_game()

# Calculate the average number of throws
average_throws = total_throws_recorded / num_simulations

# The "final equation" is the calculation of the average from the simulation data.
# We print all numbers involved in this equation as requested.
print("Results from the simulation:")
print(f"Total throws across all games: {total_throws_recorded}")
print(f"Total number of games simulated: {num_simulations}")
print("\nFinal Equation for the average:")
print(f"{total_throws_recorded} / {num_simulations} = {average_throws}")

print("\nConclusion:")
print(f"The simulated average is approximately {average_throws:.6f}.")
print(f"This value converges to the mathematical constant e, which is approximately {math.e:.6f}.")