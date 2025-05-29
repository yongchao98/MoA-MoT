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
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("IKPCINBJEUUMF")

# Round 1
cards.append('I')
cards.remove('U')
result1 = determine_result(cards, 'I')

# Round 2
cards.append('I')
cards.remove('P')
result2 = determine_result(cards, 'I')

# Round 3
cards.append('I')
cards.remove('I')
result3 = determine_result(cards, 'I')

# Round 4
cards.append('M')
cards.remove('N')
result4 = determine_result(cards, 'M')

# Round 5
cards.append('N')
result5 = determine_result(cards, 'N')

print(result5)