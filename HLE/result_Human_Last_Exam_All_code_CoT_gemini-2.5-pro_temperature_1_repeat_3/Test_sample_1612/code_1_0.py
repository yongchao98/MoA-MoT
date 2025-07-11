# A list of all batting Triple Crown wins in the American League and National League.
# Note: A player's name appears for each time they won.
al_wins = [
    "Nap Lajoie",
    "Ty Cobb",
    "Jimmie Foxx",
    "Lou Gehrig",
    "Ted Williams",
    "Ted Williams", # Won a second time
    "Mickey Mantle",
    "Frank Robinson",
    "Carl Yastrzemski",
    "Miguel Cabrera"
]

nl_wins = [
    "Paul Hines",
    "Tip O'Neill",
    "Hugh Duffy",
    "Heinie Zimmerman",
    "Rogers Hornsby",
    "Rogers Hornsby", # Won a second time
    "Chuck Klein",
    "Joe Medwick"
]

# Get the number of unique winners from each league by converting lists to sets
num_unique_al_winners = len(set(al_wins))
num_unique_nl_winners = len(set(nl_wins))

# Calculate the total number of unique winners
total_unique_winners = num_unique_al_winners + num_unique_nl_winners

# Print the final equation
print(f"Number of unique AL winners: {num_unique_al_winners}")
print(f"Number of unique NL winners: {num_unique_nl_winners}")
print(f"Total Triple Crown Winners: {num_unique_al_winners} + {num_unique_nl_winners} = {total_unique_winners}")