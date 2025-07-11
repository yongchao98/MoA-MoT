# Number of unique players who have won the batting Triple Crown
num_batting_winners = 16

# Number of unique players who have won the pitching Triple Crown
num_pitching_winners = 25

# Calculate the total number of unique Triple Crown winners
total_winners = num_batting_winners + num_pitching_winners

# Print the breakdown and the final equation
print(f"There are {num_batting_winners} unique batting Triple Crown winners.")
print(f"There are {num_pitching_winners} unique pitching Triple Crown winners.")
print(f"The total number of Triple Crown winners is the sum of these two groups.")
print(f"{num_batting_winners} + {num_pitching_winners} = {total_winners}")