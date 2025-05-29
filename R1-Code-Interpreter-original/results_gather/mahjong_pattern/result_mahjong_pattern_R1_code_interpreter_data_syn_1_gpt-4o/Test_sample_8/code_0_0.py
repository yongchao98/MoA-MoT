def check_peng(cards, new_card):
    # Check for "Peng"
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("TYHKXRZWCDIFD")

# Round 1
new_card = 'Y'
cards.append(new_card)
cards.remove('I')
if check_peng(cards, new_card):
    result_round_1 = 1
elif check_chi(cards, new_card):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
new_card = 'K'
cards.append(new_card)
cards.remove('D')
if check_peng(cards, new_card):
    result_round_2 = 1
elif check_chi(cards, new_card):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
new_card = 'I'
cards.append(new_card)
if check_peng(cards, new_card):
    result_round_3 = 1
elif check_chi(cards, new_card):
    result_round_3 = 2
else:
    result_round_3 = 0

# Output the result of the final round
print(result_round_3)