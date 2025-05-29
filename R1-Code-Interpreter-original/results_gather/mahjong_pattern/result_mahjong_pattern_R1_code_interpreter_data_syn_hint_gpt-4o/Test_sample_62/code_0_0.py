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
cards = list("TSZVJSKHGBNPL")

# Round 1
cards.append('S')
result1 = determine_result(cards, 'S')
cards.remove('S')

# Round 2
cards.append('L')
result2 = determine_result(cards, 'L')
cards.remove('N')

# Round 3
cards.append('S')
result3 = determine_result(cards, 'S')
cards.remove('L')

# Round 4
cards.append('U')
result4 = determine_result(cards, 'U')
cards.remove('S')

# Round 5
cards.append('S')
result5 = determine_result(cards, 'S')

print(result5)