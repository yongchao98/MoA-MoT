def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("QEGPWNYOAKBZN")

# Round 1
cards.append('A')
cards.remove('O')
if check_peng(cards, 'A'):
    result = 1
elif check_chi(cards, 'A'):
    result = 2
else:
    result = 0

# Round 2
cards.append('C')
cards.remove('N')
if check_peng(cards, 'C'):
    result = 1
elif check_chi(cards, 'C'):
    result = 2
else:
    result = 0

# Round 3
cards.append('R')
cards.remove('G')
if check_peng(cards, 'R'):
    result = 1
elif check_chi(cards, 'R'):
    result = 2
else:
    result = 0

# Round 4
cards.append('W')
if check_peng(cards, 'W'):
    result = 1
elif check_chi(cards, 'W'):
    result = 2
else:
    result = 0

print(result)