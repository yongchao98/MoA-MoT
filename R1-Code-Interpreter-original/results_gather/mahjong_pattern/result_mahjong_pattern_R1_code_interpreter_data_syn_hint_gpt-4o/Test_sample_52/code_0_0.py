def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for round_info in rounds:
        add_card, remove_card = round_info
        
        # Add the new card
        cards.append(add_card)
        
        # Check for "Peng"
        if cards.count(add_card) >= 3:
            result = 1
        else:
            # Check for "Chi"
            sorted_cards = sorted(cards)
            found_chi = False
            for i in range(len(sorted_cards) - 2):
                if (ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and
                    ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2):
                    found_chi = True
                    break
            result = 2 if found_chi else 0
        
        # Remove the specified card
        cards.remove(remove_card)
    
    return result

initial_cards = "HFAZRDQNCLTBA"
rounds = [('A', 'A'), ('C', 'A'), ('X', 'B'), ('U', 'Z'), ('E', None)]
final_result = determine_result(initial_cards, rounds)
print(final_result)