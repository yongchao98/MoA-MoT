def determine_result(cards, new_card):
    # Check for "Peng"
    if cards.count(new_card) >= 2:
        return 1  # Peng
    
    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if (ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and 
            ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2):
            return 2  # Chi
    
    return 0  # Pass

# Initial cards
cards = list("JVTGCAUIFBCHC")

# Round 4: Add a V card
new_card = 'V'
result = determine_result(cards, new_card)

print(result)