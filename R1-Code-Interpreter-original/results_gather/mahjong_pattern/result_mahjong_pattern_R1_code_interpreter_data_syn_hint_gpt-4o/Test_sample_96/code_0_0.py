def determine_result(cards, new_card):
    # Check for "Peng"
    if cards.count(new_card) >= 2:
        return 1  # Peng
    
    # Check for "Chi"
    sorted_cards = sorted(cards + new_card)
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi
    
    return 0  # Pass

# Initial cards
cards = list("YKUBPMMOOXAUF")

# Round 1
cards.append('L')
cards.remove('P')
result1 = determine_result(cards, 'L')

# Round 2
cards.append('M')
cards.remove('A')
result2 = determine_result(cards, 'M')

# Round 3
cards.append('O')
cards.remove('U')
result3 = determine_result(cards, 'O')

# Round 4
cards.append('I')
cards.remove('O')
result4 = determine_result(cards, 'I')

# Round 5
cards.append('M')
result5 = determine_result(cards, 'M')

print(result5)