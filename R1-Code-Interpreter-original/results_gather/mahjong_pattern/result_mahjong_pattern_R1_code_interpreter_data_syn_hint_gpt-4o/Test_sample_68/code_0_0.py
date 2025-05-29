def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("WOYZRBFFXVTQQ")

# Round 1
cards.append('X')
cards.remove('F')
if check_peng(cards, 'X'):
    result = 1
elif check_chi(cards, 'X'):
    result = 2
else:
    result = 0

# Round 2
cards.append('X')
cards.remove('V')
if check_peng(cards, 'X'):
    result = 1
elif check_chi(cards, 'X'):
    result = 2
else:
    result = 0

# Round 3
cards.append('A')
cards.remove('Y')
if check_peng(cards, 'A'):
    result = 1
elif check_chi(cards, 'A'):
    result = 2
else:
    result = 0

# Round 4
cards.append('N')
cards.remove('Z')
if check_peng(cards, 'N'):
    result = 1
elif check_chi(cards, 'N'):
    result = 2
else:
    result = 0

# Round 5
cards.append('P')
if check_peng(cards, 'P'):
    result = 1
elif check_chi(cards, 'P'):
    result = 2
else:
    result = 0

print(result)