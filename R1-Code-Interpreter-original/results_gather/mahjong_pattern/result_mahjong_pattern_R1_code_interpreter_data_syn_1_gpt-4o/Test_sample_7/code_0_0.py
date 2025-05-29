# Initial cards
cards = list("NNTMOFBJDEGAO")

# Round 1: Add 'N', Remove 'T'
cards.append('N')
cards.remove('T')

# Check for "Peng" in Round 1
if cards.count('N') >= 3:
    result_round_1 = 1
else:
    result_round_1 = 0

# Round 2: Add 'O'
cards.append('O')

# Check for "Peng" in Round 2
if cards.count('O') >= 3:
    result_round_2 = 1
else:
    result_round_2 = 0

# The result of the final round is the result of Round 2
print(result_round_2)