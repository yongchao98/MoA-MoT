# Number of unique players who have won the batting Triple Crown.
# This includes 9 in the American League and 5 in the National League.
batting_winners = 14

# Number of unique players who have won the pitching Triple Crown.
# This includes 11 in the American League and 10 in the National League.
pitching_winners = 21

# Calculate the total number of unique winners.
# No player has won both a batting and a pitching Triple Crown.
total_winners = batting_winners + pitching_winners

print(f"To find the total number of MLB Triple Crown winners, we add the unique batting winners and the unique pitching winners:")
print(f"{batting_winners} + {pitching_winners} = {total_winners}")