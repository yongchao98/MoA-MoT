def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("PQAKVDAONPQSU")

# Rounds
rounds = [
    ('A', 'Q'),
    ('F', 'A'),
    ('Q', 'V'),
    ('D', 'Q'),
    ('L', None)  # No card to remove in the last round
]

# Process each round
for new_card, remove_card in rounds:
    cards.append(new_card)
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)

# Output the result of the final round
print(result)