def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("POORAAXVJGKFE")

# Rounds
rounds = [
    ('O', 'X'),
    ('N', 'J'),
    ('Q', 'G'),
    ('A', 'O'),
    ('A', None)  # No card to remove in the last round
]

result = 0
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

print(result)