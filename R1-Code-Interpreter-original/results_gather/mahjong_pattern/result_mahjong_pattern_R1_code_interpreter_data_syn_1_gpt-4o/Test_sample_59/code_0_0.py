def check_peng(cards, new_card):
    # Check for "Peng" condition
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi" condition
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
result1 = determine_result(cards, 'F')
cards.append('F')

# Round 2
cards.remove('Q')
result2 = determine_result(cards, 'G')
cards.append('G')

# Round 3
cards.remove('R')
result3 = determine_result(cards, 'U')
cards.append('U')

# Round 4
cards.remove('E')
result4 = determine_result(cards, 'G')
cards.append('G')

# Round 5
result5 = determine_result(cards, 'F')

# Print the result of the final round
print(result5)