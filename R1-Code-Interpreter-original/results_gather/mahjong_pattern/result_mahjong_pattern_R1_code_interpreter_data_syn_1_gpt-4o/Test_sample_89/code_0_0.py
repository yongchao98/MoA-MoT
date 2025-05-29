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
cards = list("VVWOKGGWEHSYF")

# Round 1
cards.append('V')
result = determine_result(cards, 'V')
cards.remove('G')

# Round 2
cards.append('H')
result = determine_result(cards, 'H')
cards.remove('V')

# Round 3
cards.append('Q')
result = determine_result(cards, 'Q')
cards.remove('S')

# Round 4
cards.append('F')
result = determine_result(cards, 'F')
cards.remove('G')

# Round 5
cards.append('W')
result = determine_result(cards, 'W')

print(result)