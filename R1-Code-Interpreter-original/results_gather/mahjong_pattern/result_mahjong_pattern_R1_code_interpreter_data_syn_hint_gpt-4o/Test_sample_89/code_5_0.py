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
cards = list("VVWOKGGWEHSYF")

# Round 1
cards.append('V')
cards.remove('G')
result1 = determine_result(cards, 'V')

# Round 2
cards.append('H')
cards.remove('V')
result2 = determine_result(cards, 'H')

# Round 3
cards.append('Q')
cards.remove('S')
result3 = determine_result(cards, 'Q')

# Round 4
cards.append('F')
cards.remove('G')
result4 = determine_result(cards, 'F')

# Round 5
cards.append('W')
result5 = determine_result(cards, 'W')

print(result5)