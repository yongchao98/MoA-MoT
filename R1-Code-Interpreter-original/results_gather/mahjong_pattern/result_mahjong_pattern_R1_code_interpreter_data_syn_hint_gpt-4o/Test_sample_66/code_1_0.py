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
            if (index > 0 and index < len(sorted_cards) - 1 and
                ord(sorted_cards[index - 1]) == ord(new_card) - 1 and
                ord(sorted_cards[index + 1]) == ord(new_card) + 1):
                result = 2
            else:
                result = 0
    
    return result

initial_cards = "BJMGVENDQLIYI"
rounds = [('V', 'I'), ('N', 'N'), ('N', 'M'), ('N', 'L'), ('K', None)]
result = determine_result(initial_cards, rounds)
print(result)