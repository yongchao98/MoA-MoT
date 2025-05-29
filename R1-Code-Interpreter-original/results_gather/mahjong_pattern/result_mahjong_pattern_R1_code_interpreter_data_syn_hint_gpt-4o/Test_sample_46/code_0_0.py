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
cards = list("QOBSYIAVCJNYK")

# Round 1
cards.append('P')
cards.remove('C')
if check_peng(cards, 'P'):
    result = 1
elif check_chi(cards, 'P'):
    result = 2
else:
    result = 0

# Round 2
cards.append('U')
cards.remove('J')
if check_peng(cards, 'U'):
    result = 1
elif check_chi(cards, 'U'):
    result = 2
else:
    result = 0

# Round 3
cards.append('W')
cards.remove('S')
if check_peng(cards, 'W'):
    result = 1
elif check_chi(cards, 'W'):
    result = 2
else:
    result = 0

# Round 4
cards.append('Z')
cards.remove('Y')
if check_peng(cards, 'Z'):
    result = 1
elif check_chi(cards, 'Z'):
    result = 2
else:
    result = 0

# Round 5
cards.append('[')
if check_peng(cards, '['):
    result = 1
elif check_chi(cards, '['):
    result = 2
else:
    result = 0

print(result)