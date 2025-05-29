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
cards = list("QLWYBEZMPWJYF")

# Round 1
cards.append('G')
cards.remove('B')
if check_peng(cards, 'G'):
    result1 = 1
elif check_chi(cards, 'G'):
    result1 = 2
else:
    result1 = 0

# Round 2
cards.append('O')
cards.remove('J')
if check_peng(cards, 'O'):
    result2 = 1
elif check_chi(cards, 'O'):
    result2 = 2
else:
    result2 = 0

# Round 3
cards.append('Y')
cards.remove('L')
if check_peng(cards, 'Y'):
    result3 = 1
elif check_chi(cards, 'Y'):
    result3 = 2
else:
    result3 = 0

# Round 4
cards.append('W')
cards.remove('P')
if check_peng(cards, 'W'):
    result4 = 1
elif check_chi(cards, 'W'):
    result4 = 2
else:
    result4 = 0

# Round 5
cards.append('I')
if check_peng(cards, 'I'):
    result5 = 1
elif check_chi(cards, 'I'):
    result5 = 2
else:
    result5 = 0

print(result5)