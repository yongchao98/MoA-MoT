# In MLB, there are two types of Triple Crowns: one for batters and one for pitchers.
# We need to count the unique winners for each category and then sum them up.

# Number of unique players who have won the batting Triple Crown.
# This includes players from the National League, American League, and the American Association.
# Players like Rogers Hornsby and Ted Williams won it multiple times but are counted once.
num_batting_winners = 16

# Number of unique players who have won the pitching Triple Crown.
# Players like Walter Johnson, Sandy Koufax, and Roger Clemens won it multiple times but are counted once.
num_pitching_winners = 23

# Calculate the total number of unique Triple Crown winners.
total_winners = num_batting_winners + num_pitching_winners

# Print the final equation showing how the total was calculated.
print(f"Number of unique batting Triple Crown winners: {num_batting_winners}")
print(f"Number of unique pitching Triple Crown winners: {num_pitching_winners}")
print(f"Total number of unique Triple Crown winners: {num_batting_winners} + {num_pitching_winners} = {total_winners}")