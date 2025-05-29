# Initial cards
cards = list("ASCQPGYIKNVBD")

# Round 1: Add 'I', Remove 'D'
cards.append('I')
cards.remove('D')

# Check for "Peng" in Round 1
if cards.count('I') >= 3:
    round_1_result = 1
else:
    round_1_result = 0

# Round 2: Add 'V'
cards.append('V')

# Check for "Peng" in Round 2
if cards.count('V') >= 3:
    round_2_result = 1
else:
    round_2_result = 0

# Output the result of the final round
print(round_2_result)