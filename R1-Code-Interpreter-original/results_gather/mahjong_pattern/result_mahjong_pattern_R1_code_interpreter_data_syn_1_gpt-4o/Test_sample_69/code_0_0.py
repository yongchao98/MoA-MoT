def determine_result(cards, rounds):
    for add_card, remove_card in rounds:
        # Add the new card
        cards.append(add_card)
        # Remove the specified card
        cards.remove(remove_card)
        
        # Check for "Peng"
        for card in set(cards):
            if cards.count(card) >= 3:
                return 1  # "Peng"
        
        # Check for "Chi"
        sorted_cards = sorted(cards)
        for i in range(len(sorted_cards) - 2):
            if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
                return 2  # "Chi"
    
    return 0  # "Pass"

# Initial cards
initial_cards = list("QPVKTEYYUKQDC")

# Rounds: (add_card, remove_card)
rounds = [
    ('R', 'U'),
    ('Q', 'R'),
    ('Y', 'C'),
    ('K', 'K'),
    ('Q', None)  # In the last round, we don't remove a card
]

# Determine the result of the final round
result = determine_result(initial_cards, rounds)
print(result)