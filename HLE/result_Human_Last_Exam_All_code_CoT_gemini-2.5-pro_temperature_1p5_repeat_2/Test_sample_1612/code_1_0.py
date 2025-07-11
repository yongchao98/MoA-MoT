# The Triple Crown can be won by either a batter or a pitcher.
# We need to count the unique players for each category and add them together.

# Number of unique players who have won the batting Triple Crown.
# This includes players from both the American League and National League.
# Players like Ted Williams and Rogers Hornsby who won multiple times are only counted once.
num_batting_winners = 16

# Number of unique players who have won the pitching Triple Crown.
# This includes players from both leagues.
# Players like Walter Johnson, Sandy Koufax, and Lefty Grove are only counted once.
num_pitching_winners = 19

# Calculate the total number of unique Triple Crown winners.
total_winners = num_batting_winners + num_pitching_winners

# Print the final equation and the result.
print(f"Total Triple Crown Winners = Batting Winners + Pitching Winners")
print(f"Total Triple Crown Winners = {num_batting_winners} + {num_pitching_winners}")
print(f"Total = {total_winners}")