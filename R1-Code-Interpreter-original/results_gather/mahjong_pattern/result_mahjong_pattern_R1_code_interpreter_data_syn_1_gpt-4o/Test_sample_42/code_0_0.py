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
cards = list("BJUVFLMLCGKWP")

# Round 1
cards.remove('B')
cards.append('L')
result1 = determine_result(cards, 'L')

# Round 2
cards.remove('W')
cards.append('H')
result2 = determine_result(cards, 'H')

# Round 3
cards.remove('H')
cards.append('Z')
result3 = determine_result(cards, 'Z')

# Round 4
cards.remove('F')
cards.append('X')
result4 = determine_result(cards, 'X')

# Round 5
cards.append('Y')
result5 = determine_result(cards, 'Y')

print(result5)