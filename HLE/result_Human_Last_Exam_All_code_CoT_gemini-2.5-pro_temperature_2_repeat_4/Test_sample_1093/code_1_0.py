import random
import math

# Plan: We will simulate Gauss's game a large number of times to find the
# average number of throws required for the sum to exceed 1. This value
# is known to converge to Euler's number, e.

# Number of times to run the simulation for a precise average
num_simulations = 2000000

# This will store the cumulative sum of throws from all games
total_throws_recorded = 0

print(f"Running {num_simulations:,} simulations of the game...")

# Loop through each simulation
for _ in range(num_simulations):
    # Reset the sum and the throw count for the new game
    current_sum = 0.0
    num_throws = 0
    
    # Keep adding random numbers until the sum is greater than 1
    while current_sum <= 1:
        # A random number between 0 and 1, as per the problem description
        random_value = random.random()
        current_sum += random_value
        num_throws += 1
        
    # Add the number of throws from this game to our total
    total_throws_recorded += num_throws

# Calculate the average by dividing the total throws by the number of games
average_result = total_throws_recorded / num_simulations

# To meet the requirement of outputting each number in the final equation,
# we will display the components of the average calculation.
print("\n--- Calculation of the Average ---")
print("The final equation for the simulated average is:")
print(f"Total Throws / Number of Games = Average")
print(f"{total_throws_recorded} / {num_simulations} = {average_result:.7f}")

print("\n--- Conclusion ---")
print(f"The simulated average converges to: {average_result:.7f}")
print(f"This is very close to the theoretical value, Euler's number (e), which is approximately {math.e:.7f}.")
