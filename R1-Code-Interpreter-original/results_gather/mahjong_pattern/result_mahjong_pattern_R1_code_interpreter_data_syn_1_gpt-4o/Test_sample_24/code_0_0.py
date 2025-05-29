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
cards = list("TZXLGFUSYUMIK")

# Round 1
cards.append('U')
result_1 = determine_result(cards, 'U')
cards.remove('U')

# Round 2
cards.append('H')
result_2 = determine_result(cards, 'H')
cards.remove('Z')

# Round 3
cards.append('U')
result_3 = determine_result(cards, 'U')
cards.remove('X')

# Round 4
cards.append('U')
result_4 = determine_result(cards, 'U')
cards.remove('S')

# Round 5
cards.append('E')
result_5 = determine_result(cards, 'E')

print(result_5)