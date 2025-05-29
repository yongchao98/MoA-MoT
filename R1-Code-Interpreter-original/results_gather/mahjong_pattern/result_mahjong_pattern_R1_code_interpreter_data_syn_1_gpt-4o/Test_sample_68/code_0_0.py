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
cards = list("WOYZRBFFXVTQQ")

# Round 1
cards.append('X')
cards.remove('F')
result1 = determine_result(cards, 'X')

# Round 2
cards.append('X')
cards.remove('V')
result2 = determine_result(cards, 'X')

# Round 3
cards.append('A')
cards.remove('Y')
result3 = determine_result(cards, 'A')

# Round 4
cards.append('N')
cards.remove('Z')
result4 = determine_result(cards, 'N')

# Round 5
cards.append('P')
result5 = determine_result(cards, 'P')

print(result5)