# Initial cards
cards = list("IKTNTHZGBSFRF")

def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if adding the new card forms a sequence with any two other cards
    cards.append(new_card)
    sorted_cards = sorted(set(cards))
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Round 1
new_card = 'F'
cards.append(new_card)
cards.remove('S')
if check_peng(cards, new_card):
    result_round_1 = 1
elif check_chi(cards, new_card):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
new_card = 'F'
cards.append(new_card)
cards.remove('F')
if check_peng(cards, new_card):
    result_round_2 = 1
elif check_chi(cards, new_card):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
new_card = 'F'
cards.append(new_card)
cards.remove('F')
if check_peng(cards, new_card):
    result_round_3 = 1
elif check_chi(cards, new_card):
    result_round_3 = 2
else:
    result_round_3 = 0

# Round 4
new_card = 'T'
cards.append(new_card)
cards.remove('H')
if check_peng(cards, new_card):
    result_round_4 = 1
elif check_chi(cards, new_card):
    result_round_4 = 2
else:
    result_round_4 = 0

# Round 5
new_card = 'A'
cards.append(new_card)
if check_peng(cards, new_card):
    result_round_5 = 1
elif check_chi(cards, new_card):
    result_round_5 = 2
else:
    result_round_5 = 0

print(result_round_5)