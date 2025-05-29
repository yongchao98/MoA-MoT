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
cards = list("XDMTLKZWUJWHR")

# Round 1
cards.append('W')
result1 = determine_result(cards, 'W')
cards.remove('L')

# Round 2
cards.append('Y')
result2 = determine_result(cards, 'Y')
cards.remove('W')

# Round 3
cards.append('P')
result3 = determine_result(cards, 'P')
cards.remove('J')

# Round 4
cards.append('X')
result4 = determine_result(cards, 'X')
cards.remove('W')

# Round 5
cards.append('V')
result5 = determine_result(cards, 'V')

# Output the result of the final round
print(result5)