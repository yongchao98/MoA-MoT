import collections

# A list of all MLB batting Triple Crown wins.
# Each tuple contains the year, player's name, and league.
triple_crown_wins = [
    (1878, "Paul Hines", "NL"),
    (1887, "Tip O'Neill", "AA"),
    (1901, "Nap Lajoie", "AL"),
    (1909, "Ty Cobb", "AL"),
    (1912, "Heinie Zimmerman", "NL"),
    (1922, "Rogers Hornsby", "NL"),
    (1925, "Rogers Hornsby", "NL"),
    (1933, "Jimmie Foxx", "AL"),
    (1933, "Chuck Klein", "NL"),
    (1934, "Lou Gehrig", "AL"),
    (1937, "Joe Medwick", "NL"),
    (1942, "Ted Williams", "AL"),
    (1947, "Ted Williams", "AL"),
    (1956, "Mickey Mantle", "AL"),
    (1966, "Frank Robinson", "AL"),
    (1967, "Carl Yastrzemski", "AL"),
    (2012, "Miguel Cabrera", "AL"),
]

# Use a set to get the unique player names to avoid double-counting.
unique_winners = sorted(list(set(player for year, player, league in triple_crown_wins)))

print("The MLB Triple Crown winners are:")
for winner in unique_winners:
    print(winner)

# Calculate the total number of unique winners.
num_winners = len(unique_winners)

# Create and print the equation as requested.
equation_parts = ["1"] * num_winners
equation_str = " + ".join(equation_parts)

print(f"\nCounting each unique winner:")
print(f"{equation_str} = {num_winners}")

print(f"\nIn total, there are {num_winners} unique Triple Crown winners in MLB history.")
