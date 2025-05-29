def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for add_card, remove_card in rounds:
        # Add the new card
        cards.append(add_card)
        
        # Check for "Peng"
        if cards.count(add_card) >= 3:
            return 1  # Peng
        
        # Check for "Chi"
        sorted_cards = sorted(cards)
        for i in range(len(sorted_cards) - 2):
            if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
                return 2  # Chi
        
        # Remove the specified card
        cards.remove(remove_card)
    
    return 0  # Pass

initial_cards = "SDCBJJBLPELSN"
rounds = [('S', 'S'), ('L', 'C'), ('S', 'B'), ('H', 'S'), ('S', None)]
result = determine_result(initial_cards, rounds)
print(result)