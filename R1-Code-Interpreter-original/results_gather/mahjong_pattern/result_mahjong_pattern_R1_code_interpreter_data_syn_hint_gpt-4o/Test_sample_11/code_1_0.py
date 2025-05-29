# Initial cards
cards = list("TWWZFUUVIMPYB")

def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    cards_set = set(cards)
    if (chr(ord(new_card) - 1) in cards_set and chr(ord(new_card) - 2) in cards_set) or \
       (chr(ord(new_card) + 1) in cards_set and chr(ord(new_card) + 2) in cards_set) or \
       (chr(ord(new_card) - 1) in cards_set and chr(ord(new_card) + 1) in cards_set):
        return True
    return False

# Round 1
cards.append('T')
cards.remove('V')
if check_peng(cards, 'T'):
    result_round_1 = 1
elif check_chi(cards, 'T'):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
cards.append('R')
cards.remove('W')
if check_peng(cards, 'R'):
    result_round_2 = 1
elif check_chi(cards, 'R'):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
cards.append('I')
if check_peng(cards, 'I'):
    result_round_3 = 1
elif check_chi(cards, 'I'):
    result_round_3 = 2
else:
    result_round_3 = 0

# Output the result of the final round
print(result_round_3)