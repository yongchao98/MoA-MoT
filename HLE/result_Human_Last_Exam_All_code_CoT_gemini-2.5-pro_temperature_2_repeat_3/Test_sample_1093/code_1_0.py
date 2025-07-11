import random

# The number of times Gauss's game is simulated.
# A higher number leads to a more accurate average.
num_simulations = 1000000

# This variable will store the sum of the recorded numbers from all games.
total_throws_recorded = 0

# Run the simulation num_simulations times.
for _ in range(num_simulations):
    
    # Reset for a new game.
    sum_of_measurements = 0.0
    num_throws = 0
    
    # Steps 1-3: Keep adding measurements until the sum is larger than 1.
    while sum_of_measurements <= 1:
        # Generate a random number between 0 and 1.
        measurement = random.random()
        sum_of_measurements += measurement
        num_throws += 1
        
    # Add the recorded number of throws for this game to the total.
    total_throws_recorded += num_throws

# Calculate the final average.
average = total_throws_recorded / num_simulations

# As requested, output each number in the final equation for the average.
print(f"The average is calculated as follows:")
print(f"Total Throws Recorded / Number of Simulations = Average")
print(total_throws_recorded, "/", num_simulations, "=", average)