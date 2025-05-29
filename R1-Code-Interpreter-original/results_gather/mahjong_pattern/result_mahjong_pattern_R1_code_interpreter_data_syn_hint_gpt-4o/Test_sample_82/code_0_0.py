def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1
    elif check_chi(cards, new_card):
        return 2
    else:
        return 0

# Initial cards
cards = list("RESMSABPXGYTF")

# Round 1
cards.append('G')
cards.remove('P')
result1 = determine_result(cards, 'G')

# Round 2
cards.append('T')
cards.remove('R')
result2 = determine_result(cards, 'T')

# Round 3
cards.append('S')
cards.remove('A')
result3 = determine_result(cards, 'S')

# Round 4
cards.append('H')
cards.remove('S')
result4 = determine_result(cards, 'H')

# Round 5
cards.append('B')
result5 = determine_result(cards, 'B')

print(result5)