def determine_result(cards, new_card):
    # Add the new card
    cards += new_card
    
    # Check for "Peng"
    for card in set(cards):
        if cards.count(card) >= 3:
            return 1  # Peng
    
    # Check for "Chi"
    sorted_cards = sorted(cards)
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi
    
    return 0  # Pass

# Initial cards
cards = "WGLRQQYUWVLXZ"

# Round 1
cards = cards.replace('W', '', 1)  # Remove W
result1 = determine_result(cards, 'Q')

# Round 2
cards = cards.replace('V', '', 1)  # Remove V
result2 = determine_result(cards, 'Q')

# Round 3
cards = cards.replace('W', '', 1)  # Remove W
result3 = determine_result(cards, 'Z')

# Round 4
cards = cards.replace('X', '', 1)  # Remove X
result4 = determine_result(cards, 'Q')

# Round 5
result5 = determine_result(cards, 'Q')

print(result5)