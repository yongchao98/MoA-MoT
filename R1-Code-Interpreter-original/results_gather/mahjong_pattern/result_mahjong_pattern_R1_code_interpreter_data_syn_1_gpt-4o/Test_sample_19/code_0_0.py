def determine_result(cards, new_card):
    # Check for "Peng"
    if cards.count(new_card) >= 2:
        return 1  # Peng

    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("JVTGCAUISFBCR")

# Round 1
cards.append('C')
result1 = determine_result(cards, 'C')
cards.remove('C')

# Round 2
cards.append('H')
result2 = determine_result(cards, 'H')
cards.remove('S')

# Round 3
cards.append('C')
result3 = determine_result(cards, 'C')
cards.remove('R')

# Round 4
cards.append('V')
result4 = determine_result(cards, 'V')

print(result4)