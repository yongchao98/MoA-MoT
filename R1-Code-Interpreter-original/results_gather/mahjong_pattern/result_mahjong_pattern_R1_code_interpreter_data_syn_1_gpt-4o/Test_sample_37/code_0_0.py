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
cards = list("GASNDQYSZLRCW")

# Round 1
cards.remove('W')
cards.append('H')
result1 = determine_result(cards, 'H')

# Round 2
cards.remove('N')
cards.append('E')
result2 = determine_result(cards, 'E')

# Round 3
cards.remove('Z')
cards.append('I')
result3 = determine_result(cards, 'I')

# Round 4
cards.remove('D')
cards.append('P')
result4 = determine_result(cards, 'P')

# Round 5
cards.append('I')
result5 = determine_result(cards, 'I')

print(result5)