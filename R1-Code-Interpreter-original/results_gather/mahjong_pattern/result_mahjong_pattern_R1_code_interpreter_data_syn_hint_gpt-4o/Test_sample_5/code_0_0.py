# Initial set of cards
cards = list("YQILCUCZEAXNK")

# Round 1: Add 'H', remove 'E'
cards.append('H')
cards.remove('E')

# Check for "Peng" in Round 1
round_1_result = 0  # Default to "Pass"
if cards.count('H') >= 3:
    round_1_result = 1  # "Peng"

# Round 2: Add 'C'
cards.append('C')

# Check for "Peng" in Round 2
round_2_result = 0  # Default to "Pass"
if cards.count('C') >= 3:
    round_2_result = 1  # "Peng"

# Output the result of the final round
print(round_2_result)