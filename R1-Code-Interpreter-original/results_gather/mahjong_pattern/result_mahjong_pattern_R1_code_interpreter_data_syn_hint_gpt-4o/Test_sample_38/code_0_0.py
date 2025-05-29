def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for add_card, remove_card in rounds:
        cards.append(add_card)
        cards.remove(remove_card)
        
        # Check for "Peng"
        if cards.count(add_card) >= 3:
            result = 1
        else:
            # Check for "Chi"
            sorted_cards = sorted(set(cards))
            for i in range(len(sorted_cards) - 2):
                if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
                    result = 2
                    break
            else:
                result = 0
    
    return result

initial_cards = "VQTWUJOAYYJGE"
rounds = [('J', 'A'), ('X', 'V'), ('Y', 'J'), ('Z', 'Y'), ('Z', '')]

final_result = determine_result(initial_cards, rounds)
print(final_result)