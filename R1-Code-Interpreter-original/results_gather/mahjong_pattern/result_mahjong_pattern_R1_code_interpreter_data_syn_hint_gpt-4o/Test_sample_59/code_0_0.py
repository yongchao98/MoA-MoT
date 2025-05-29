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
cards = list("DTRKLKQUNXOME")

# Round 1
cards.remove('K')
cards.append('F')
result1 = determine_result(cards, 'F')

# Round 2
cards.remove('Q')
cards.append('G')
result2 = determine_result(cards, 'G')

# Round 3
cards.remove('R')
cards.append('U')
result3 = determine_result(cards, 'U')

# Round 4
cards.remove('E')
cards.append('G')
result4 = determine_result(cards, 'G')

# Round 5
result5 = determine_result(cards, 'F')

print(result5)