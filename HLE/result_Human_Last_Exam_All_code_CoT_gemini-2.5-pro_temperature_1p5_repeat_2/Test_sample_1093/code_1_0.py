import random
import math

# Set the number of times we repeat the experiment.
# A higher number will give a result closer to the theoretical value.
num_trials = 1000000

total_count = 0

# Loop through each trial.
for _ in range(num_trials):
    sum_of_measurements = 0.0
    count = 0
    # Keep adding measurements until the sum is larger than 1.
    while sum_of_measurements <= 1:
        # Generate a random number between 0 and 1, simulating the needle throw.
        measurement = random.random()
        sum_of_measurements += measurement
        count += 1
    
    # Record the number of throws for this trial.
    total_count += count

# Calculate the average over all trials.
average = total_count / num_trials

print(f"Simulating the game for {num_trials} trials.")
print("The average number of throws required for the sum to exceed 1 is calculated below.")
# The final equation with each number printed out
print(f"{total_count} / {num_trials} = {average}")
print(f"\nThis value is an approximation for Euler's number, e.")
print(f"The actual value of e is {math.e}")
