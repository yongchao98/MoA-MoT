def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for add_card, remove_card in rounds:
        cards.append(add_card)
        if remove_card:
            cards.remove(remove_card)
        
        # Check for "Peng"
        for card in set(cards):
            if cards.count(card) >= 3:
                return 1  # Peng
        
        # Check for "Chi"
        unique_cards = sorted(set(cards))
        for i in range(len(unique_cards) - 2):
            if (ord(unique_cards[i+1]) == ord(unique_cards[i]) + 1 and
                ord(unique_cards[i+2]) == ord(unique_cards[i]) + 2):
                return 2  # Chi
    
    return 0  # Pass

initial_cards = "QNWSHPNZFKSUF"
rounds = [('F', 'N'), ('D', 'S'), ('R', None)]
result = determine_result(initial_cards, rounds)
print(result)