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
cards = list("FSNGLAWBLONTQ")

# Round 1
cards.remove('A')
cards.append('H')

# Round 2
cards.remove('B')
cards.append('Y')

# Round 3
cards.remove('Q')
cards.append('K')

# Round 4
cards.remove('S')
cards.append('L')

# Round 5
cards.append('N')

# Determine the result for the final round
result = determine_result(cards, 'N')
print(result)