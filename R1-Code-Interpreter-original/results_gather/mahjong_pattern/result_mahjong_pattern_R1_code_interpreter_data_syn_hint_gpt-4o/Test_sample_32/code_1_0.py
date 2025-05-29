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
        if check_peng(cards, new_card):
            result = 1
        elif check_chi(cards, new_card):
            result = 2
        else:
            result = 0
        if remove_card is not None and remove_card in cards:
            cards.remove(remove_card)
    return result

initial_cards = "HAHPSQTMLCDXB"
rounds = [('D', 'D'), ('L', 'Q'), ('W', 'A'), ('D', 'C'), ('L', None)]
result = play_game(initial_cards, rounds)
print(result)