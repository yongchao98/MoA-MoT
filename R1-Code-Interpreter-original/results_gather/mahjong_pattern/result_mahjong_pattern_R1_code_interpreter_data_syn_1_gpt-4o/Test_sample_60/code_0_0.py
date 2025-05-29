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
cards = list("GAWRMZMPOWJHV")

# Round 1
cards.append('J')
cards.remove('R')
result1 = determine_result(cards, 'J')

# Round 2
cards.append('W')
cards.remove('Z')
result2 = determine_result(cards, 'W')

# Round 3
cards.append('P')
cards.remove('P')
result3 = determine_result(cards, 'P')

# Round 4
cards.append('M')
cards.remove('O')
result4 = determine_result(cards, 'M')

# Round 5
cards.append('W')
result5 = determine_result(cards, 'W')

print(result5)