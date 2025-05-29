# Initial cards
cards = list("IKTNTHZGBSFRF")

# Round 1
cards.append('F')
cards.remove('S')
if cards.count('F') >= 3:
    result_round_1 = 1
else:
    result_round_1 = 0

# Round 2
cards.append('F')
cards.remove('F')
if cards.count('F') >= 3:
    result_round_2 = 1
else:
    result_round_2 = 0

# Round 3
cards.append('F')
cards.remove('F')
if cards.count('F') >= 3:
    result_round_3 = 1
else:
    result_round_3 = 0

# Round 4
cards.append('T')
cards.remove('H')
if cards.count('T') >= 3:
    result_round_4 = 1
else:
    # Check for Chi
    sorted_cards = sorted(set(cards))
    chi_found = False
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            chi_found = True
            break
    result_round_4 = 2 if chi_found else 0

# Round 5
cards.append('A')
if cards.count('A') >= 3:
    result_round_5 = 1
else:
    # Check for Chi
    sorted_cards = sorted(set(cards))
    chi_found = False
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            chi_found = True
            break
    result_round_5 = 2 if chi_found else 0

print(result_round_5)