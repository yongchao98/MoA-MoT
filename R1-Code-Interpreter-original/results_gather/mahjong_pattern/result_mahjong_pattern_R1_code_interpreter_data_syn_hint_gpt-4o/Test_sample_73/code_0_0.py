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
cards = list("QNGKSTZTDHCOE")

# Round 1
cards.remove('Q')
cards.append('P')
result1 = determine_result(cards, 'P')

# Round 2
cards.remove('G')
cards.append('T')
result2 = determine_result(cards, 'T')

# Round 3
cards.remove('K')
cards.append('E')
result3 = determine_result(cards, 'E')

# Round 4
cards.remove('T')
cards.append('R')
result4 = determine_result(cards, 'R')

# Round 5
cards.append('U')
result5 = determine_result(cards, 'U')

print(result5)