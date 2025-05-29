# Initial cards
cards = list("VQFIUPNXEKOEJ")

# Round 1
cards.append('E')  # Add E
cards.remove('V')  # Remove V
# Check for "Peng"
round_1_result = 1 if cards.count('E') >= 3 else 0

# Round 2
cards.append('P')  # Add P
# Check for "Peng"
round_2_result = 1 if cards.count('P') >= 3 else 0

# The final result is the result of the last round
print(round_2_result)