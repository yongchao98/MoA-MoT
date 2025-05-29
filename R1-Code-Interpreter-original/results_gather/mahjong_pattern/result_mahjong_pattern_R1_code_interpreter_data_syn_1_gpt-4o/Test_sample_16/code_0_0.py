def check_peng(cards, new_card):
    # Check if there are two identical cards in hand and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("LBHKYWOMRFDSI")

# Round 1
cards.append('N')
cards.remove('B')
if check_peng(cards, 'N'):
    result_round_1 = 1
elif check_chi(cards, 'N'):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
cards.append('V')
cards.remove('D')
if check_peng(cards, 'V'):
    result_round_2 = 1
elif check_chi(cards, 'V'):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
cards.append('N')
cards.remove('F')
if check_peng(cards, 'N'):
    result_round_3 = 1
elif check_chi(cards, 'N'):
    result_round_3 = 2
else:
    result_round_3 = 0

# Round 4
cards.append('Q')
if check_peng(cards, 'Q'):
    result_round_4 = 1
elif check_chi(cards, 'Q'):
    result_round_4 = 2
else:
    result_round_4 = 0

# Output the result of the final round
print(result_round_4)