import random

# Set the number of games to simulate. A higher number will give a result
# closer to the theoretical value.
num_simulations = 2000000

total_recorded_numbers = 0

# Loop through each game simulation.
for _ in range(num_simulations):
    # Step 4: Reset sum for the new game.
    current_sum = 0.0
    num_throws = 0
    
    # Step 1-3: Keep adding random numbers until the sum is greater than 1,
    # and count how many numbers were needed.
    while current_sum <= 1:
        # Step 1-2: Generate a random number between 0 and 1.
        measurement = random.random()
        current_sum += measurement
        num_throws += 1
    
    # Add the recorded count for this game to the total.
    total_recorded_numbers += num_throws

# Calculate the average over all simulated games.
average = total_recorded_numbers / num_simulations

# Print the final equation with each number.
# The average should be very close to e (~2.71828).
print("The average is derived from the following calculation:")
print(f"{total_recorded_numbers} / {num_simulations} = {average}")