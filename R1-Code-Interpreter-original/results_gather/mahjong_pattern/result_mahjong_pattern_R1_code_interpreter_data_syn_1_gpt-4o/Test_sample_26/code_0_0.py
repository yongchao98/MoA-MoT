def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + new_card)
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
cards = list("NDMMGOARIIQEZ")

# Round 1
cards.remove('I')
result1 = determine_result(cards, 'I')
cards.append('I')

# Round 2
cards.remove('M')
result2 = determine_result(cards, 'I')
cards.append('I')

# Round 3
cards.remove('D')
result3 = determine_result(cards, 'P')
cards.append('P')

# Round 4
cards.remove('A')
result4 = determine_result(cards, 'O')
cards.append('O')

# Round 5
result5 = determine_result(cards, 'S')

print(result5)