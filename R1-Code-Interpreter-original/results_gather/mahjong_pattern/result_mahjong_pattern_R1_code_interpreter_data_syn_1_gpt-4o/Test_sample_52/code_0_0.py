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
cards = list("HFAZRDQNCLTBA")

# Round 1
cards.append('A')
result1 = determine_result(cards, 'A')
cards.remove('A')

# Round 2
cards.append('C')
result2 = determine_result(cards, 'C')
cards.remove('A')

# Round 3
cards.append('X')
result3 = determine_result(cards, 'X')
cards.remove('B')

# Round 4
cards.append('U')
result4 = determine_result(cards, 'U')
cards.remove('Z')

# Round 5
cards.append('E')
result5 = determine_result(cards, 'E')

print(result5)