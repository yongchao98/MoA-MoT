def check_peng(cards, new_card):
    # Check for "Peng" condition
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi" condition
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("RNMKZVUXGWFZD")

# Round 1
cards.append('B')
cards.remove('K')
if check_peng(cards, 'B'):
    result_round_1 = 1
elif check_chi(cards, 'B'):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
cards.append('W')
cards.remove('U')
if check_peng(cards, 'W'):
    result_round_2 = 1
elif check_chi(cards, 'W'):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
cards.append('Y')
if check_peng(cards, 'Y'):
    result_round_3 = 1
elif check_chi(cards, 'Y'):
    result_round_3 = 2
else:
    result_round_3 = 0

# Output the result of the final round
print(result_round_3)