# Initial cards
cards = list("TWWZFUUVIMPYB")

# Round 1
cards.append('T')
cards.remove('V')
# Check for Peng
if cards.count('T') >= 3:
    result_round_1 = 1
else:
    # Check for Chi
    sorted_cards = sorted(set(cards))
    chi_found = False
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            chi_found = True
            break
    result_round_1 = 2 if chi_found else 0

# Round 2
cards.append('R')
cards.remove('W')
# Check for Peng
if cards.count('R') >= 3:
    result_round_2 = 1
else:
    # Check for Chi
    sorted_cards = sorted(set(cards))
    chi_found = False
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            chi_found = True
            break
    result_round_2 = 2 if chi_found else 0

# Round 3
cards.append('I')
# Check for Peng
if cards.count('I') >= 3:
    result_round_3 = 1
else:
    # Check for Chi
    sorted_cards = sorted(set(cards))
    chi_found = False
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            chi_found = True
            break
    result_round_3 = 2 if chi_found else 0

# Output the result of the final round
print(result_round_3)