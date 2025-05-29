def check_peng(cards, new_card):
    # Check for "Peng"
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi"
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
cards = list("GRCPONYJVADUA")

# Round 1
cards.append('A')
cards.remove('C')
result1 = determine_result(cards, 'A')

# Round 2
cards.append('N')
cards.remove('Y')
result2 = determine_result(cards, 'N')

# Round 3
cards.append('W')
cards.remove('A')
result3 = determine_result(cards, 'W')

# Round 4
cards.append('Q')
cards.remove('D')
result4 = determine_result(cards, 'Q')

# Round 5
cards.append('S')
result5 = determine_result(cards, 'S')

# Final result of the last round
print(result5)