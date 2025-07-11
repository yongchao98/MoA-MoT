import random
import math

# --- Introduction ---
print("J.C.F. Gauss's game is a classic probability problem.")
print("The average number of random throws (from 0 to 1) needed for their sum to exceed 1 converges to a famous mathematical constant: Euler's number, e.")
print("\n--- Part 1: The Theoretical Value (e) ---")
print("The value of e can be calculated from the infinite series:")
print("e = 1/0! + 1/1! + 1/2! + 1/3! + ...\n")
print("Let's calculate this sum for the first 15 terms to see it converge:")

e_approximation = 0.0
# Fulfilling the request to "output each number in the final equation"
# by showing the terms of the series for e.
for n in range(15):
    term = 1.0 / math.factorial(n)
    e_approximation += term
    print(f"  Term n={n:2d}: 1/{n}! = {term:.8f}")

print("-" * 45)
print(f"Sum of first 15 terms: {e_approximation:.8f}")
print(f"Actual value of e:     {math.e:.8f}")
print("=" * 45)

# --- Part 2: Simulation of the Game ---
print("\n--- Part 2: Simulating the Game ---")
print("Now, we'll run a simulation to verify this theoretical result.")

# Number of times to run the entire experiment
num_simulations = 2000000

# This will store the sum of the results (number of throws) from all simulations
total_throws_recorded = 0

# Loop to run the simulation many times
for i in range(num_simulations):
    sum_of_measurements = 0.0
    num_throws = 0
    # Keep throwing the needle until the sum is greater than 1
    while sum_of_measurements <= 1:
        # A single throw gives a random number between 0 and 1
        random_value = random.random()
        sum_of_measurements += random_value
        num_throws += 1
    
    # Record the number of throws it took for this game
    total_throws_recorded += num_throws

# Calculate the average number of throws over all simulations
average_result = total_throws_recorded / num_simulations

print(f"\nSimulation finished after running the game {num_simulations:,} times.")
print("\nThe experimental average is calculated by dividing the total number of throws by the number of games played.")

# Fulfilling the request to "output each number in the final equation"
# by showing the final average calculation.
print("\nFinal Equation from Simulation:")
print(f"{total_throws_recorded} (Total Throws) / {num_simulations} (Games) = {average_result:.8f}")
<<<e>>>