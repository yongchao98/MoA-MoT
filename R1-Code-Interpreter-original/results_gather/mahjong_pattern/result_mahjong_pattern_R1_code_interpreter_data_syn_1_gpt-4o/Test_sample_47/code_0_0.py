def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    card_values = sorted([ord(c) for c in cards] + [ord(new_card)])
    for i in range(len(card_values) - 2):
        if card_values[i+1] == card_values[i] + 1 and card_values[i+2] == card_values[i] + 2:
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
cards = list("XQFTKUZIXJGBD")

# Round 1
cards.remove('Z')
cards.append('L')
result1 = determine_result(cards, 'L')

# Round 2
cards.remove('F')
cards.append('D')
result2 = determine_result(cards, 'D')

# Round 3
cards.remove('K')
cards.append('W')
result3 = determine_result(cards, 'W')

# Round 4
cards.remove('L')
cards.append('X')
result4 = determine_result(cards, 'X')

# Round 5
result5 = determine_result(cards, 'I')

print(result5)