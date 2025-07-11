import random

# Set the number of times to simulate the game.
# A higher number will produce an average closer to the theoretical value.
num_simulations = 2000000

# This list will store the count from each game.
results = []

print(f"Simulating Gauss's game {num_simulations} times...")
print("This may take a moment.")

# Loop to run the simulation the specified number of times.
for _ in range(num_simulations):
    
    # Step 1-3: Sum random numbers until the total is > 1 and record the count.
    current_sum = 0.0
    count = 0
    while current_sum <= 1:
        # Generate a random number between 0 and 1.
        random_number = random.random()
        current_sum += random_number
        count += 1
    
    # Add the count for this game to our list of results.
    results.append(count)

# Calculate the total sum of counts and the number of games played.
total_sum_of_counts = sum(results)
number_of_games = len(results)

# Calculate the average.
average_result = total_sum_of_counts / number_of_games

print("\nSimulation complete.")
print("The average of the recorded numbers is the sum of all counts divided by the number of games.")
print("The final equation is:")
print(f"{total_sum_of_counts} / {number_of_games} = {average_result}")
print("\nThis value is an approximation of Euler's number, e (≈ 2.71828).")
