import random

# The number of times to run the simulation to get a stable average.
# A higher number will give a result closer to the theoretical value.
num_simulations = 2000000

total_throws_recorded = 0

# Step 4: Reset the sum and repeat the game many times.
for _ in range(num_simulations):
    
    # This block represents one full game (Steps 1-3)
    sum_of_measurements = 0.0
    number_of_throws = 0
    
    # Step 1-2: Add measurements until the sum is larger than 1.
    while sum_of_measurements <= 1:
        # Get a random number between 0 and 1.
        measurement = random.uniform(0, 1)
        sum_of_measurements += measurement
        number_of_throws += 1
    
    # Step 3: Record the number of throws for this game.
    total_throws_recorded += number_of_throws

# Calculate the average of the recorded numbers.
# The average is the total number of throws divided by the number of games played.
average_result = total_throws_recorded / num_simulations

# The "final equation" is the calculation of the average from the simulation data.
# Here we output each number used in that final calculation.
print("Final Equation from Simulation:")
print(f"Average = Total Throws / Number of Simulations")
print(f"Average = {total_throws_recorded} / {num_simulations}")
print(f"The simulated average is: {average_result}")
print("\nAs shown by the simulation, this value converges to the theoretical value, e (Euler's number).")
print(f"e ≈ 2.71828")
