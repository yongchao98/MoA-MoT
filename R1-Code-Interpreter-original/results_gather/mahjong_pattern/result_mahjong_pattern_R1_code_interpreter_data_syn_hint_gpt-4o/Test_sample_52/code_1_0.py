def determine_result(initial_cards, rounds):
    from collections import Counter
    
    # Initialize card counts
    card_counts = Counter(initial_cards)
    
    for add_card, remove_card in rounds:
        # Add the new card
        card_counts[add_card] += 1
        
        # Check for "Peng"
        if card_counts[add_card] >= 3:
            result = 1
        else:
            # Check for "Chi"
            unique_cards = sorted(card_counts.keys())
            found_chi = False
            for i in range(len(unique_cards) - 2):
                if (ord(unique_cards[i+1]) == ord(unique_cards[i]) + 1 and
                    ord(unique_cards[i+2]) == ord(unique_cards[i]) + 2):
                    found_chi = True
                    break
            result = 2 if found_chi else 0
        
        # Remove the specified card
        if remove_card:
            card_counts[remove_card] -= 1
            if card_counts[remove_card] == 0:
                del card_counts[remove_card]
    
    return result

initial_cards = "HFAZRDQNCLTBA"
rounds = [('A', 'A'), ('C', 'A'), ('X', 'B'), ('U', 'Z'), ('E', None)]
final_result = determine_result(initial_cards, rounds)
print(final_result)