def determine_result(initial_cards, rounds):
    cards = list(initial_cards)
    
    for round_info in rounds:
        new_card, remove_card = round_info
        cards.append(new_card)
        
        if remove_card is not None:
            if remove_card in cards:
                cards.remove(remove_card)
        
        # Check for "Peng"
        if cards.count(new_card) >= 3:
            result = 1
        else:
            # Check for "Chi"
            sorted_cards = sorted(set(cards))
            index = sorted_cards.index(new_card)
            chi_found = False
            for i in range(len(sorted_cards) - 2):
                if (ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and
                    ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2):
                    if new_card in sorted_cards[i:i+3]:
                        chi_found = True
                        break
            if chi_found:
                result = 2
            else:
                result = 0
    
    return result

initial_cards = "BJMGVENDQLIYI"
rounds = [('V', 'I'), ('N', 'N'), ('N', 'M'), ('N', 'L'), ('K', None)]
result = determine_result(initial_cards, rounds)
print(result)