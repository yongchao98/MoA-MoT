def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for add_card, remove_card in rounds:
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

initial_cards = "HAHPSQTMLCDXB"
rounds = [('D', 'D'), ('L', 'Q'), ('W', 'A'), ('D', 'C'), ('L', None)]

# For the last round, we don't remove any card, so we handle it separately
final_result = determine_result(initial_cards, rounds[:-1])
# Add the last card without removing any
final_result = determine_result(initial_cards + rounds[-1][0], [])

print(final_result)