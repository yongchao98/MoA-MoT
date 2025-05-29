def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    cards = sorted(cards + [new_card])
    for i in range(len(cards) - 2):
        if ord(cards[i+1]) == ord(cards[i]) + 1 and ord(cards[i+2]) == ord(cards[i]) + 2:
            return True
    return False

def game_result(initial_cards, rounds):
    cards = list(initial_cards)
    for i, (new_card, remove_card) in enumerate(rounds):
        if check_peng(cards, new_card):
            result = 1
        elif check_chi(cards, new_card):
            result = 2
        else:
            result = 0
        cards.append(new_card)
        if remove_card in cards:
            cards.remove(remove_card)
    return result

initial_cards = "TZXLGFUSYUMIK"
rounds = [('U', 'U'), ('H', 'Z'), ('U', 'X'), ('U', 'S'), ('E', None)]

# For the last round, we only add a card and do not remove any
final_result = game_result(initial_cards, rounds)
print(final_result)