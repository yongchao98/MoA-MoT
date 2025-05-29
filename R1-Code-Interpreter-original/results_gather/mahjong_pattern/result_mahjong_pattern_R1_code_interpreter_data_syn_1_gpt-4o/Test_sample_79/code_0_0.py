def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
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
cards = list("OQVJEDCGCGBPO")

# Rounds
rounds = [
    ('T', 'B'),
    ('R', 'E'),
    ('O', 'R'),
    ('S', 'C'),
    ('G', None)  # No card to remove in the final round
]

# Process each round
for new_card, remove_card in rounds:
    result = determine_result(cards, new_card)
    if remove_card:
        cards.remove(remove_card)
    cards.append(new_card)

# The result of the final round
print(result)