def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("PCHWEHUWNAGXO")

# Round 1
cards.append('P')
cards.remove('W')
if check_peng(cards, 'P'):
    result = 1
elif check_chi(cards, 'P'):
    result = 2
else:
    result = 0

# Round 2
cards.append('Z')
cards.remove('N')
if check_peng(cards, 'Z'):
    result = 1
elif check_chi(cards, 'Z'):
    result = 2
else:
    result = 0

# Round 3
cards.append('I')
cards.remove('H')
if check_peng(cards, 'I'):
    result = 1
elif check_chi(cards, 'I'):
    result = 2
else:
    result = 0

# Round 4
cards.append('H')
cards.remove('H')
if check_peng(cards, 'H'):
    result = 1
elif check_chi(cards, 'H'):
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