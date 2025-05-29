def determine_result(cards, new_card):
    # Check for "Peng"
    for card in set(cards):
        if cards.count(card) == 2 and new_card == card:
            return 1  # Peng

    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("SWDMIGCRTNLQV")

# Round 1
cards.append('T')
cards.remove('I')
result1 = determine_result(cards, 'T')

# Round 2
cards.append('O')
cards.remove('N')
result2 = determine_result(cards, 'O')

# Round 3
cards.append('T')
cards.remove('T')
result3 = determine_result(cards, 'T')

# Round 4
cards.append('T')
cards.remove('M')
result4 = determine_result(cards, 'T')

# Round 5
cards.append('X')
result5 = determine_result(cards, 'X')

print(result5)