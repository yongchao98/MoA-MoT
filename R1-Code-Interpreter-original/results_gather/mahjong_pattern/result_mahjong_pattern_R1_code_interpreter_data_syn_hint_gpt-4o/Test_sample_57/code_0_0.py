def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

def play_game(initial_cards, rounds):
    cards = list(initial_cards)
    for new_card, remove_card in rounds:
        cards.append(new_card)
        if remove_card:
            cards.remove(remove_card)
        
        if check_peng(cards, new_card):
            result = 1
        elif check_chi(cards, new_card):
            result = 2
        else:
            result = 0
    
    return result

initial_cards = "PQAKVDAONPQSU"
rounds = [('A', 'Q'), ('F', 'A'), ('Q', 'V'), ('D', 'Q'), ('L', None)]

final_result = play_game(initial_cards, rounds)
print(final_result)