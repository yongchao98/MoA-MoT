def determine_result(cards, new_card):
    # Check for Peng
    if cards.count(new_card) >= 2:
        return 1  # Peng
    
    # Check for Chi
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi
    
    return 0  # Pass

# Initial cards
cards = list("XMPGJPLYNIJBH")

# Round 1
cards.remove('H')
cards.append('F')
result1 = determine_result(cards, 'F')

# Round 2
cards.remove('I')
cards.append('Z')
result2 = determine_result(cards, 'Z')

# Round 3
cards.remove('F')
cards.append('R')
result3 = determine_result(cards, 'R')

# Round 4
cards.remove('Z')
cards.append('X')
result4 = determine_result(cards, 'X')

# Round 5
cards.append('J')
result5 = determine_result(cards, 'J')

print(result5)