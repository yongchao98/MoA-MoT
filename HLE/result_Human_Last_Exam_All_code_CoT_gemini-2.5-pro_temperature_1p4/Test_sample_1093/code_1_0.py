import random

def run_one_game():
    """
    Simulates one round of the game as described by J.C.F. Gauss.
    It generates random numbers between 0 and 1 and adds them up
    until the sum is larger than 1.

    Returns:
        int: The number of random numbers (throws) it took for the sum to exceed 1.
    """
    sum_of_measurements = 0.0
    number_of_throws = 0
    while sum_of_measurements <= 1:
        # Step 1-2: Generate a random number between 0 and 1
        measurement = random.random()
        sum_of_measurements += measurement
        number_of_throws += 1
    # Step 3: Record the number
    return number_of_throws

# --- Main Simulation ---

# Step 4: Reset and repeat a large number of times to find the average
num_simulations = 2000000
total_throws_across_simulations = 0

for i in range(num_simulations):
    total_throws_across_simulations += run_one_game()

# Calculate the average
average_result = total_throws_across_simulations / num_simulations

# The problem states to output each number in the final equation.
# Here is the final calculation based on the simulation.
print("--- Simulation Results ---")
print(f"Total number of throws recorded: {total_throws_across_simulations}")
print(f"Total number of games simulated: {num_simulations}")
print("Final average calculation:")
print(f"{total_throws_across_simulations} / {num_simulations} = {average_result}")
print("\nThe simulation shows that the average of the recorded numbers converges to this value.")

# The theoretical value this simulation approaches is Euler's number, e.
# import math
# print(f"\nFor reference, the exact theoretical value is e ≈ {math.e}")
