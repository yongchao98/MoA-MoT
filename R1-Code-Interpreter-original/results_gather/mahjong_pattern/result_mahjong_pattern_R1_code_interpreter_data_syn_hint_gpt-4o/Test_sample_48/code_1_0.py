def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for round_info in rounds:
        new_card, remove_card = round_info
        cards.append(new_card)
        
        if remove_card is not None and remove_card in cards:
            cards.remove(remove_card)
        
        # Check for "Peng"
        if cards.count(new_card) >= 3:
            result = 1
        else:
            # Check for "Chi"
            sorted_cards = sorted(cards)
            chi_found = False
            for i in range(len(sorted_cards) - 2):
                if (ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and
                    ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2):
                    chi_found = True
                    break
            result = 2 if chi_found else 0
    
    return result

initial_cards = "WBZGBEZNTIWFQ"
rounds = [('B', 'W'), ('N', 'E'), ('U', 'Z'), ('X', 'T'), ('Y', None)]
final_result = determine_result(initial_cards, rounds)
print(final_result)