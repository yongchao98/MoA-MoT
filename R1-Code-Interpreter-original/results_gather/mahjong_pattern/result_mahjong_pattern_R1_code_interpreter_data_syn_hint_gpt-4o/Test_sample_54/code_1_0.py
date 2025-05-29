# Initial cards
cards = list("IKPCINBJEUUMF")

# Function to check for "Peng"
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Function to check for "Chi"
def check_chi(cards, new_card):
    cards_with_new = cards + [new_card]
    cards_with_new.sort()
    for i in range(len(cards_with_new) - 2):
        if (ord(cards_with_new[i+1]) == ord(cards_with_new[i]) + 1 and
            ord(cards_with_new[i+2]) == ord(cards_with_new[i]) + 2):
            return True
    return False

# Round 1
cards.append('I')
cards.remove('U')
if check_peng(cards, 'I'):
    result = 1
elif check_chi(cards, 'I'):
    result = 2
else:
    result = 0

# Round 2
cards.append('I')
cards.remove('P')
if check_peng(cards, 'I'):
    result = 1
elif check_chi(cards, 'I'):
    result = 2
else:
    result = 0

# Round 3
cards.append('I')
cards.remove('I')
if check_peng(cards, 'I'):
    result = 1
elif check_chi(cards, 'I'):
    result = 2
else:
    result = 0

# Round 4
cards.append('M')
cards.remove('N')
if check_peng(cards, 'M'):
    result = 1
elif check_chi(cards, 'M'):
    result = 2
else:
    result = 0

# Round 5
cards.append('N')
if check_peng(cards, 'N'):
    result = 1
elif check_chi(cards, 'N'):
    result = 2
else:
    result = 0

print(result)