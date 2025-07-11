import random
import math

def run_simulation(num_trials):
    """
    Simulates Gauss's game for a given number of trials to find the
    average number of throws needed for the sum of random numbers to exceed 1.
    """
    total_throws_across_all_trials = 0

    # Loop for the number of specified trials
    for _ in range(num_trials):
        current_sum = 0.0
        throws_for_this_trial = 0
        # Keep throwing until the sum is larger than 1
        while current_sum <= 1:
            # Generate a random number between 0 and 1
            random_value = random.random()
            current_sum += random_value
            throws_for_this_trial += 1
        
        # Add the number of throws from this trial to the total
        total_throws_across_all_trials += throws_for_this_trial

    # Calculate the average
    average_throws = total_throws_across_all_trials / num_trials
    return average_throws

# --- Main Execution ---

# Set a large number of trials for a more accurate result
num_trials = 2000000

# Run the simulation
simulation_average = run_simulation(num_trials)

# --- Explanation and Final Answer ---

print(f"Running a simulation of the game for {num_trials:,} trials...")
print(f"The average number of throws from the simulation is: {simulation_average:.6f}\n")

print("This simulation demonstrates a classic mathematical result.")
print("The average number of throws converges to Euler's number, 'e'.")
print("The theoretical value is derived from the sum of the probabilities for the number of throws needed.")
print("\nThe final equation for the exact average is the Taylor series for e:")
print("Average = 1/0! + 1/1! + 1/2! + 1/3! + 1/4! + ...")
print("\nHere are the numbers that make up the start of this infinite sum:")

# Calculate and print each number in the final equation
term_0 = 1 / math.factorial(0)
term_1 = 1 / math.factorial(1)
term_2 = 1 / math.factorial(2)
term_3 = 1 / math.factorial(3)
term_4 = 1 / math.factorial(4)
term_5 = 1 / math.factorial(5)
partial_sum = term_0 + term_1 + term_2 + term_3 + term_4 + term_5

print(f"Term 1 (1/0!): {term_0}")
print(f"Term 2 (1/1!): {term_1}")
print(f"Term 3 (1/2!): {term_2}")
print(f"Term 4 (1/3!): {term_3:.6f}")
print(f"Term 5 (1/4!): {term_4:.6f}")
print(f"Term 6 (1/5!): {term_5:.6f}")
print(f"The sum continues indefinitely, and its sum is e ≈ {math.e:.6f}")
print(f"The sum of just these first 6 terms is: {partial_sum:.6f}")

<<<e>>>